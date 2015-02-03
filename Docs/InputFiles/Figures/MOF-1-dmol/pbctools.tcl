#####################################################################################
#                                                                                   #
# PBCTOOLS:                                                                         #
#                                                                                   #
# Various tools to set up and display periodic boundary conditions.                 #
# See documentation in the script headers.                                          #
#                                                                                   #
# Author:                                                                           #
#   Jan Saam                                                                        #
#   Institute of Biochemistry, Charite                                              #
#   Monbijoustr. 2, Berlin                                                          #
#   Germany                                                                         #
#   saam@charite.de                                                                 #
#                                                                                   #
# Thanks to Cameron Mura <cmura@mccammon.ucsd.edu> for valuable additions.          #       
#                                                                                   #
#####################################################################################
package provide pbctools 1.2

###########################################################
# Set the periodic cell for all frames dimensions given   #
# in VMD format.                                          #
###########################################################
proc set_unitcell_vmd {a b c {alpha 90.0} {beta 90.0} {gamma 90.0} {molid top}} {
    set numframes [molinfo $molid get numframes]

    for {set i 0} {$i < $numframes} {incr i} {
        molinfo $molid set frame $i
        molinfo $molid set a $a
        molinfo $molid set b $b
        molinfo $molid set c $c
        molinfo $molid set alpha $alpha
        molinfo $molid set beta $beta
        molinfo $molid set gamma $gamma
    }
}

###########################################################
# Get the periodic cell dimensions given in VMD format.   #
###########################################################
proc get_unitcell_vmd {{molid top}} {

   puts "a = [molinfo $molid get a]"
   puts "b = [molinfo $molid get b]"
   puts "c = [molinfo $molid get c]"
   puts "alpha = [molinfo $molid get alpha]  (angle b-c)"
   puts "beta  = [molinfo $molid get beta ]  (angle a-c)"
   puts "gamma = [molinfo $molid get gamma]  (angle a-b, x-y plane)"
   return [molinfo $molid get {a b c alpha beta gamma}]
}

 
#####################################################
# Get the periodic cell dimensions.                 #
#####################################################

proc get_unitcell { {molid top} } {
   # In molinfo the length of the cell vectors and the angles between 
   # them are saved. $a is assumed to point in x-direction.
   # Notes: a, b, c are side lengths of the unit cell
   # alpha = angle between b and c
   # beta  = angle between a and c
   # gamma = angle between a and b
   
   set cell [molinfo $molid get {a b c alpha beta gamma}]

   return [unitcell_vmd2namd $cell]
}



###############################################################
# set_unitcell                                                #
#                                                             #
# Set the periodic cell dimensions, given the PBC vectors     #
# such as they appear in the XST/XSC or in NAMD config files  #
# Vector A should parallel to x-axis, otherwise VMD           #
# cannot display periodic cells properly!                     #
#                                                             #
# Options:                                                    #
#    -alignx    rotate unitcell vector A into the x-axis.     #
#               The system will also be rotated accordingly.  #
#    -draw      draw the PBC-vectors.                         #
#                                                             #
###############################################################

proc set_unitcell { cell args } {
   # Defaults
   set molid  top
   set alignx 0;   # Automaticly align A with the x-axis (0/1)?
   set draw   0;   # Draw the cell (0/1)?
   set retval 0;   # Return value will be 1 if the coordinates were rotated

   # Scan for single options
   set argnum 0
   set arglist $args
   foreach i $args {
      if {$i=="-alignx"}  then {
         set alignx 1
         set arglist [lreplace $arglist $argnum $argnum]
	 continue
      }
      if {$i=="-draw"}  then {
         set draw 1
         set arglist [lreplace $arglist $argnum $argnum]
	 continue
      }
      incr argnum
   }

   # Scan for options with one argument
   foreach {i j} $arglist {
      if {$i=="-molid"}   { set molid  $j }
   }

   set A [lindex $cell 0]
   set B [lindex $cell 1]
   set C [lindex $cell 2]
   set rot {}

   # In molinfo the length of the cell vectors and the angles between 
   # them are saved. $a is assumed to point in x-direction.
   if {[expr abs([vecdot [vecnorm $A] {1 0 0}]-1.0)]>0.000001} {
      if {$alignx} {
	 # Compute transformation matrix to rotate A into x-axis
	 set rot [transvecinv $A]
	 set retval 1
      } else {
	 append msg "Vector A is not parallel to x-axis, \n"
	 append msg "VMD cannot display periodic cells properly!\n"
	 append msg "Use option '-alignx' to rotate unitcell vector A into\n"
	 append msg "the x-axis. The system will also be rotated accordingly.\n"
	 append msg "You can specify '-xaxis auto' to choose it automatically."
	 error $msg
      }
   }

   if {$draw} { draw_pbc_vectors $cell }

   if {[llength $rot]} {
      set A [coordtrans $rot $A]
      set B [coordtrans $rot $B]
      set C [coordtrans $rot $C]
      if {$draw} { draw_pbc_vectors [list $A $B $C {0 0 0}] }
      set sel [atomselect $molid all]
      $sel move $rot
   }

   set a [veclength $A]
   set b [veclength $B]
   set c [veclength $C]
   set alpha 1
   set beta  1
   set gamma 1

   # Note: Between VMD 1.8.2 and 1.8.3 the definition of the unitcell
   # parameters changed which is why we have to check the version
   if {[string compare "1.8.3" [vmdinfo version]]>0} {
      # puts "VMD version <= 1.8.2"
      set gamma [vecangle $B $C]
      set beta  [vecangle $A $C]
      set alpha [vecangle $A $B]
   } else {
      # puts "VMD version > 1.8.2 (including 1.8.3aXX)"
      set alpha [vecangle $B $C]
      set beta  [vecangle $A $C]
      set gamma [vecangle $A $B]
   }

   #puts "set_unitcell (namd format): $cell"   
   #puts "set_unitcell (vmd format) : $a $b $c $alpha $beta $gamma"
   molinfo $molid set a $a
   molinfo $molid set b $b
   molinfo $molid set c $c
   molinfo $molid set alpha $alpha
   molinfo $molid set beta $beta
   molinfo $molid set gamma $gamma

   return $retval
}


#####################################################
# Set the periodic cell dimensions for each frame.  #
# using the info from an xst file                   #
# Options:                                          #
# -first  $f   first frame to process               #
# -last   $f   last frame to process                #
# -molid  $m   molecule ID                          #
# -stride $n   Read only every n-th line from xst   #
#              (This is useful when dcd was loaded  #
#              with a stride                        #
# -log    $l   log to write to                      #
# -noskipfirst do not skip the first line of the    #
#              XST file                             #
# -step2frame  conversion factor between step num   #
#              in XST file and frame num in DCDs    #
#              (useful when loading multiple XSTs   #
#              and want to avoid over-writing info  #
#              for earlier frames by having unique  #
#              step <-> frame mapping)              #
#####################################################

proc set_unitcell_xst { xstfile args } {
   # Defaults
   set first  0
   set last   "last"
   set stride 1
   set molid  "top"
   set skipfirst 1
   set step2frame 0
   set log    {}

   # Scan for single options
   set argnum 0
   set arglist $args
   foreach i $args {
      if {$i=="-noskipfirst"}  then {
         set skipfirst 0
         set arglist [lreplace $arglist $argnum $argnum]
	 continue
      }
      incr argnum
   }

   # Parse the arguments
   foreach {i j} $arglist {
      if {$i=="-first"}   { set first     $j }
      if {$i=="-last"}    { set last      $j }
      if {$i=="-stride"}  { set stride    $j }
      if {$i=="-molid"}   { set molid     $j }
      if {$i=="-log"}     { set log       $j }
      if {$i=="-step2frame"} { set step2frame       $j }
   }

   if { $molid=="top" }  then { set molid [ molinfo top ] }

   set numf [molinfo $molid get numframes]
   if { $last=="last" }   then {
      set last [expr $numf-1]
   }

   if {! [file exists $xstfile]} {
      error "set_unitcell_xst: Didn't find XST file $xstfile"
   }
   set fd [open "$xstfile" r]

   # Save the current frame number
   set curframe [ molinfo $molid get frame ]

   set warn    0;  # If the axis where rotated, $warn is increased
   set dt      0
   set time    0
   set numline 0
   set frame   $first
   while {![eof $fd]} {
      set line [gets $fd]
      if {[string first \# $line]==-1 && [llength $line]>0} {
         # The first line is omitted, because xst info starts at frame 0 
         # while dcd record starts at frame 1.
	 if {$skipfirst && $numline==0} { 
	    if {[llength $log]} { puts $log "Skipping first entry" }
	    incr numline
	    continue 
	 }

	 if {!($numline%$stride) && $frame<=$last} {
	    # Get the time
	    set oldtime $time;
	    set olddt $dt
	    set time [lrange $line 0 0]
	    # Get PBC vectors
	    set v1   [lrange $line 1 3]
	    set v2   [lrange $line 4 6]
	    set v3   [lrange $line 7 9]
	    set ori  [lrange $line 10 12]
	    set cell [list $v1 $v2 $v3 $ori]

	    # Check if the number of timesteps per frame changed
	    set dt [expr $time-$oldtime];
	    if {!$numline==1 && $dt!=$olddt && [llength $log]} {
	       puts $log "\nWARNING Stepsize in XST changed! dt=$dt, olddt=$olddt\n"
	    }

		# if provided, use conversion factor for times > 0:
		if {$step2frame && $time} {
			if {$stride != 1} {
				set effectiveframe [expr round( $time * $step2frame / $stride ) - 1]
			} else {
				set effectiveframe [expr $time * $step2frame - 1]
			}
		} else {
			set effectiveframe $frame
		}
		#DB puts "time = $time / effectiveframe = $effectiveframe"
	    # Proceed to next frame
	    molinfo $molid set frame $effectiveframe


	    # Now actually set the unit cell
	    incr warn [set_unitcell $cell -molid $molid -alignx]

		#DB puts "set_unitcell_xst: $frame $time $cell"

	    if {[llength $log]} {
	       puts $log "set_unitcell_xst: $frame $time $cell"
	    }
	    incr frame
	 }
	 incr numline
      }
   }

   molinfo $molid set frame $curframe
   close $fd
   if {$warn && [llength $log]} {
      puts $log "WARNING: The coordinates have been rotated in $warn frames"
      puts $log "         so that unitcell vector A is parallel to the x-axis!"
   }
}


#########################################################
# Flip alpha and gamma in the periodic cell dimensions. #
# Between VMD 1.8.2 and 1.8.3 the definition of the     #
# unitcell angles changed such that vmd 1.8.3 displays  #
# DCDs with included unitcell info wrong.               #
# If you run flip_unitcell_angles after loading the     #
# trajectory it will be displayed correctly.            #
#########################################################

proc flip_unitcell_angles { {molid top} } {
   set alpha [molinfo $molid get gamma]
   set beta  [molinfo $molid get beta]
   set gamma [molinfo $molid get alpha]
   molinfo $molid set alpha $alpha
   molinfo $molid set beta  $beta
   molinfo $molid set gamma $gamma
}


###################################################################
# Transform VMD style unit cell info into NAMD unit cell vectors  #
###################################################################

proc unitcell_vmd2namd { vmdcell args } {
   set verbose 0

   foreach i $args {
      if {$i=="-v"}  then {
         set verbose 1
         set arglist [lreplace $arglist $argnum $argnum]
	 continue
      }
      incr argnum
   }

   set a     [lindex $vmdcell 0]
   set b     [lindex $vmdcell 1]
   set c     [lindex $vmdcell 2]
   set alpha [lindex $vmdcell 3]
   set beta  [lindex $vmdcell 4]
   set gamma [lindex $vmdcell 5]

   # The following is taken from VMD Timestep.C
   # void Timestep::get_transforms(Matrix4 &a, Matrix4 &b, Matrix4 &c)

   # A will lie along the positive x axis.
   # B will lie in the x-y plane
   # The origin will be (0,0,0).

   # a, b, c are side lengths of the unit cell
   # alpha = angle between b and c
   # beta  = angle between a and c
   # gamma = angle between a and b

   set A {}; set B {}; set C {};

   # Note: Between VMD 1.8.2 and 1.8.3 the definition of the unitcell
   # parameters changed which is why we have to check the version
   if {[string compare "1.8.3" [vmdinfo version]]>0} {
      # puts "VMD version <= 1.8.2"
      set alphar [deg2rad $gamma];  # swapped!
      set betar  [deg2rad $beta];
      set gammar [deg2rad $alpha];  # swapped!

      set cosAB  [expr cos($alphar)];
      set sinAB  [expr sin($alphar)];
      set cosAC  [expr cos($betar)];
      set cosBC  [expr cos($gammar)];
      
      set Ax $a
      set Bx [expr $b*$cosAB]
      set By [expr $b*$sinAB]
      set Cx [expr $c*$cosAC]
      set Cy [expr ($b*$c*$cosBC-$Bx*$Cx)/$By]
      set Cz [expr sqrt($c*$c-$Cx*$Cx-$Cy*$Cy)]
      
      set A  "$Ax 0.0 0.0"
      set B  "$Bx $By 0.0"
      set C  "$Cx $Cy $Cz"
      
      set phi [vecangle {0 0 1} $C]
      set Cl [expr $c/cos([deg2rad $phi])]
      set C [vecscale $Cl [vecnorm $C]]
   } else {
      # puts "VMD version > 1.8.2 (including 1.8.3aXX)"
      set cosBC [expr cos([deg2rad $alpha])]
      set sinBC [expr sin([deg2rad $alpha])]
      set cosAC [expr cos([deg2rad $beta])]
      set cosAB [expr cos([deg2rad $gamma])]
      set sinAB [expr sin([deg2rad $gamma])]
      
      set Ax $a
      set Bx [expr $b*$cosAB]
      set By [expr $b*$sinAB]
      
      # If sinAB is zero, then we can't determine C uniquely since it's defined
      # in terms of the angle between A and B.
      if {$sinAB>0} {
	 set Cx [expr $cosAC]
	 set Cy [expr ($cosBC - $cosAC * $cosAB) / $sinAB]
	 set Cz [expr sqrt(1.0 - $Cx*$Cx - $Cy*$Cy)]
      } else {
	 set Cx 0.0
	 set Cy 0.0
	 set Cz 0.0
      }
      
      set A  "$Ax 0.0 0.0"
      set B  "$Bx $By 0.0"
      set C  "$Cx $Cy $Cz"
      set C [vecscale $C $c]
   }

   set cella "[precn 1 $A]"
   set cellb "[precn 1 $B]"
   set cellc "[precn 1 $C]"

   if {$verbose} {
      puts "A=$cella; a=[veclength $cella]"
      puts "B=$cellb; b=[veclength $cellb]"
      puts "C=$cellc; c=[veclength $cellc]"
      puts "alpha = [vecangle $cellb $cellc]"
      puts "beta  = [vecangle $cella $cellc]"
      puts "gamma = [vecangle $cella $cellb]"
   }

   return [list $cella $cellb $cellc {0 0 0}]
}


###########################################################
# Read unitcells from xst or xsc file and return a list   #
# of vectors [list $cell1 $cell2 $cell3 $cellorigin].     #
###########################################################

proc read_xst {file} {
   set cell {}
   set cell1 {}
   set cell2 {}
   set cell3 {}
   set cello {}

   set fd [open $file r]
      while {![eof $fd]} {
         set line [gets $fd]
         if {[string first \# $line]==-1 && [llength $line]>0} {
            # Get PBC vectors
	    set cell1 [precn 2 [lrange $line 1 3]]
	    set cell2 [precn 2 [lrange $line 4 6]]
	    set cell3 [precn 2 [lrange $line 7 9]]
	    set cello [precn 2 [lrange $line 10 12]]
	    lappend cell [list $cell1 $cell2 $cell3 $cello]
         }
      }
      close $fd
   return $cell
}



######################################################
# Writes the PBC dimensions to a xsc file.           #
######################################################

proc write_xsc { cell filename } {
   set fd [open "$filename" w]
   puts $fd "# NAMD extended system configuration file. GENERATED BY pbc.tcl"
   puts $fd "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
   set a   [lindex $cell 0]
   set b   [lindex $cell 1]
   set c   [lindex $cell 2]
   set ori [lindex $cell 3]
   set strain {0 0 0 0 0 0}; # Strain rate is is ignored.
   puts $fd "0 $a $b $c $ori $strain"
   close $fd
}


######################################################
# Computes the size of the rectangularwater box that #
# is needed to cut the given unitcell out of the     #
# system.                                            #
######################################################

proc waterboxsize { cell } {		
   set v1  [lindex $cell 0]
   set v2  [lindex $cell 1]
   set v3  [lindex $cell 2]
   set ori [lindex $cell 3]
   set v123 [vecadd $ori $v1 $v2 $v3]
   set size [vecsub $v123 $ori]
   set low  [vecsub $ori [vecscale 0.5 $size]]
   set high [vecadd $low $size]

   draw color orange
   draw box $low $high 3
   return "{$low} {$high}"
}


######################################################
# Cuts the geometry specified by $cell out of $sel   #
# and saves it in a pdb file.                        #
######################################################

proc cut_pbc_cell { sel cell filename {offset 0.0}} {
   set v1  [lindex $cell 0]
   set v2  [lindex $cell 1]
   set v3  [lindex $cell 2]
   set ori [lindex $cell 3]
   set xsize [lindex $v1 0]
   set xori  [lindex $ori 0]
   set alpha [vecangle $v1 $v2]

   # You can adjust the cutting limits by an offset:
   set offset [expr double($offset)/sin([deg2rad $alpha])]
   # slope
   set m [precn 3 [expr tan([deg2rad $alpha])]]
   # x(y=0)
   set x0 [expr $xori+0.5*$xsize]
   # intercept
   set b  [precn 3 [expr $m*($x0-$offset)]  ]

   set sel [atomselect top "([$sel text]) and not same residue as y>(x*$m+$b)"]
   set sel [atomselect top "([$sel text]) and not same residue as y<(x*$m-$b)"]
   set molid [$sel molid]
   showsel $sel top name Points

   $sel writepdb $filename
}


######################################################
# Draws PBC cell dimesions                           #
# (a box containing the molecule).                   #
# This is a shortcut for 'draw_pbc_vectors -box'     #
######################################################

proc draw_pbc_cell { args } {
   if {[llength $args]==0} {
      error "Usage: draw_pbc_cell $cell \[-color $color\]"
   }
   eval draw_pbc_vectors $args -box
}


######################################################
# Draws PBC cell vectors.                            #
# Options:                                           #
# -box      a unitcell box containing the molecule   #
#           will be drawn.                           #
# -color    color of the box                         #
######################################################

proc draw_pbc_vectors { cell args } {
    set box 0
    set positive 0
    set color yellow

    # Scan for single options
    set argnum 0
    set arglist $args
    foreach i $args {
	if {$i=="-box"}  then {
	    set box 1
	    set arglist [lreplace $arglist $argnum $argnum]
	    continue
	}
	if {$i=="-positive"}  then {
	    set positive 1
	    set arglist [lreplace $arglist $argnum $argnum]
	    continue
	}
	incr argnum
    }

    # Scan for options with argument
    foreach {i j} $arglist {
	if {$i=="-color"}    then { set color $j }
    }
    
    if {[llength $cell]!=4} {
	error "Usage: draw_pbc_vectors $cell \[-color $color\]"
    }
    
    set v1  [lindex $cell 0]
    set v2  [lindex $cell 1]
    set v3  [lindex $cell 2]
    set ori [lindex $cell 3]
    if {$box} {
	if {!$positive} then {
	    set ori [vecsub $ori [vecscale 0.5 [vecadd $v1 $v2 $v3]]]
	} else {
	    set ori { 0.0 0.0 0.0 }
	}
    }
    set p1   [vecadd $ori $v1]
    set p2   [vecadd $ori $v2]
    set p3   [vecadd $ori $v3]
    set p12  [vecadd $p1 $v2]
    set p13  [vecadd $p1 $v3]
    set p23  [vecadd $p2 $v3]
    set p123 [vecadd $p12 $v3]
    if {!$box} {
	draw color orange
	draw sphere $ori radius 1
	draw color red
	draw arrow $ori $p1 1
	draw color green
	draw arrow $ori $p2 1
	draw color blue
	draw arrow $ori $p3 1
    }
    draw color $color
    draw line $ori $p1 width 2
    draw line $ori $p2 width 2
    draw line $ori $p3 width 2
    draw line $p1 $p12 width 2
    draw line $p2 $p12 width 2
    draw line $p1 $p13 width 2
    draw line $p3 $p13 width 2
    draw line $p2 $p23 width 2
    draw line $p3 $p23 width 2
    draw line $p12 $p123 width 2 
    draw line $p13 $p123 width 2
    draw line $p23 $p123 width 2
}


#########################################################
# Draws the minmax box containing the selection $sel.   #
#########################################################

proc draw_minmax_box { sel } {
   set box [measure minmax $sel]
   set min [lindex $box 0]
   set max [lindex $box 1]
   set size [vecsub $max $min]
   set v1  "[lindex $size 0]  0  0"
   set v2  "0  [lindex $size 1]  0"
   set v3  "0  0  [lindex $size 2]"
   set p1   [vecadd $min $v1]
   set p2   [vecadd $min $v2]
   set p3   [vecadd $min $v3]
   set p12  [vecadd $p1 $v2]
   set p13  [vecadd $p1 $v3]
   set p23  [vecadd $p2 $v3]
   set p123 [vecadd $p12 $v3]
   draw line $min $p1
   draw line $min $p2
   draw line $min $p3
   draw line $p1 $p12
   draw line $p2 $p12
   draw line $p1 $p13
   draw line $p3 $p13
   draw line $p2 $p23
   draw line $p3 $p23
   draw line $p12 $max
   draw line $p13 $max
   draw line $p23 $max
}


######################################################
# Draws the PBC hexagon encompassing the molecule.   #
######################################################

proc draw_pbc_hexagon {p V} {
   set p1 [lindex $p 0]
   set p2 [lindex $p 1]
   set p3 [lindex $p 2]
   set V1 [lindex $V 0]
   set V2 [lindex $V 1]
   set V3 [lindex $V 2]

   draw color green
   draw line $p1 $p2
   draw line $p2 $p3
   draw line $p1 [vecinvert $p3]
   draw line [vecinvert $p1] $p3
   draw line [vecinvert $p2] [vecinvert $p1]
   draw line [vecinvert $p2] [vecinvert $p3]

   draw color red
   draw line {0 0 0} $V1
   draw color orange
   draw line {0 0 0} $V2
   draw color yellow
   draw line {0 0 0} $V3

   # Set visual center to origin
   trans origin {0 0 0}
   display update
}

######################################################
# Draws the PBC hexagon encompassing the molecule.   #
# Debugging version: Also draws construction lines.  #
######################################################

proc draw_pbc_hexagon_debug {p angle1 angle2 D} {
   # Works in x,y-plane only!
   set V1 [lindex $V 0]
   set V2 [lindex $V 1]
   set V3 [lindex $V 2]
   draw color red
   draw line {0 0 0} $V1
   draw color orange
   draw line {0 0 0} $V2
   draw color yellow
   draw line {0 0 0} $V3

   draw color purple
   draw line {0 0 0} [lindex $D 0]
   draw line {0 0 0} [lindex $D 1]
   draw line {0 0 0} [lindex $D 2]

   set p1 [lindex $p 0]
   set p2 [lindex $p 1]
   set p3 [lindex $p 2]

   # Length of the sides
   set s1 [veclength [vecsub $p2 $p1]]
   set s2 [veclength [vecsub $p2 $p3]]
   set s3 [expr [veclength [vecsub $p1 [vecinvert $p3]]]]

   # Geometry of the wedge that has to be subtracted from the
   # area of the parallel epiped
   set gamma [expr 0.5*3.14159265-$angle1]
   set p [expr $s2 * sin($gamma)]
   set h [expr $s2 * cos($gamma)]
   set q [expr $h/tan(-($angle2))]

   draw color yellow
   set p4 "[expr [lindex $p2 0]+[sign $s1]*($p+$q)] 0 [expr $d1*0.5]"
   set p5 "[expr [lindex $p2 0]+[sign $s1]*$p]      0 [expr $d1*0.5]"
   draw line $p2 $p4
   draw line $p3 $p5 
   draw line $p3 $p4
}


######################################################
# Draws the boundaries of the PBC cell.              #
######################################################

proc draw_pbc_boundary { p1 p2 p3 cellc } {
   set h [vecscale $cellc 0.5]
   set p4 [vecinvert $p1]
   set p5 [vecinvert $p2]
   set p6 [vecinvert $p3]
   set p1 [vecadd $p1 $h]
   set p2 [vecadd $p2 $h]
   set p3 [vecadd $p3 $h]
   set p4 [vecadd $p4 $h]
   set p5 [vecadd $p5 $h]
   set p6 [vecadd $p6 $h]
   draw color green
   draw line $p1 $p2
   draw line $p2 $p3
   draw line $p1 $p6
   draw line $p4 $p3
   draw line $p5 $p4
   draw line $p5 $p6
   set q1 [vecadd $p1 [vecinvert $cellc]]
   set q2 [vecadd $p2 [vecinvert $cellc]]
   set q3 [vecadd $p3 [vecinvert $cellc]]
   set q4 [vecadd $p4 [vecinvert $cellc]]
   set q5 [vecadd $p5 [vecinvert $cellc]]
   set q6 [vecadd $p6 [vecinvert $cellc]]
   draw line $q1 $q2
   draw line $q2 $q3
   draw line $q1 $q6
   draw line $q4 $q3
   draw line $q5 $q4
   draw line $q5 $q6

   draw line $p1 $q1
   draw line $p2 $q2
   draw line $p3 $q3
   draw line $p4 $q4
   draw line $p5 $q5
   draw line $p6 $q6
}


# Other helper functions:
# =======================

###################################################
# Computes the angle between two vectors x and y  #
###################################################

proc vecangle {x y} {
    if {[llength $x] != [llength $y]} {
        error "vecangle needs arrays of the same size: $x : $y"
    }
    if {[llength $x]==0 || [llength $y]==0} {
        error "vecangle: zero length vector: [llength $x] : [llength $y]"
    }
    # Compute scalar-produt
    set dot 0
    foreach t1 $x t2 $y {
        set dot [expr $dot + $t1 * $t2]
    }
    set rr [rad2deg [expr (acos($dot/([veclength $x] * [veclength $y])))]]

    return $rr
}


###################################################
# Returns number or list with n digit precision   #
###################################################

proc precn { n list } {
   if {[llength $list]>1} {
      set newlist {}
      foreach val $list {
         lappend newlist [precn $n $val]
      }
      return $newlist
   } else {
      return [expr pow(10, -$n)*(round($list*pow(10, $n)))]
   }
}


###################################################
# Draws an arrow.                                 #
###################################################

#proc vmd_draw_arrow {mol start end {rad 1} {res 6}} {
    # an arrow is made of a cylinder and a cone
#    set middle [vecadd $start [vecscale 0.85 [vecsub $end $start]]]
#    graphics $mol cone $middle $end radius [expr $rad*2.0] resolution $res
#    graphics $mol cylinder $start $middle radius $rad resolution $res
    #puts "$middle $end"
#}

proc deg2rad { deg } {
  return [expr ($deg/180.0*3.14159265)]
}

proc rad2deg { rad } {
  return [expr ($rad/3.14159265)*180.0]
}


# Excerpt of communication with John Stone, explaining different behavior
# of VMD versions 1.8.2 and 1.8.3:
# -----------------------------------------------------------------------

# Jan,
#   Bad news, it looks like the Charmm guys have changed the DCD
# format again.  Jim and I just spent some time sorting through various
# DCD permutations, and it would appear that when we updated the DCD plugin
# for VMD to correctly handle some Charmm formats, we didn't catch all of
# the bizarre cases.  So, the likely cause of your problem is a difference
# between what NAMD is writing out and what the DCD reader in VMD is
# interpreting, caused by the changes that have been made to retain
# Charmm compatibility....  Sigh...

# Anyway, we'll get the DCD issue sorted out soon, hopefully that's
# the real culprit here, and not the code in VMD itself.

#   John

# On Mon, Dec 20, 2004 at 05:13:22PM -0600, John Stone wrote:

# >>
# >> Jan,
# >>   Yes, the PBC code was changed in VMD 1.8.3 to correct incorrect
# >> PBC cell rendering.  The new code now follows the standard
# >> crystallographic conventions for the meaning of a/b/c/alpha/beta/gamma
# >> as published in Sands' book "Introduction to Crystallography", and by
# >> testing against quite a few datasets in several other file formats.
# >> This is the same convention used by the PDB file format.
# >>
# >> Regarding your ability to define A vectors that are not aligned
# >> with the X axis, this an obvious limitation of using the crystallographically
# >> oriented PBC description method rather than the use of cell basis vectors.
# >> I intend to update the VMD plugin interfaces for the next version, adding
# >> support for the use of cell basis vectors, which should address your
# >> concern, though not in the short term.  (I have to update the plugin interface,
# >> the Tcl/Python interfaces, and come up with a way so that scripts can determine
# >> if the cell basis can be represented as a/b/c/alpha/beta/gamma or not rather
# >> than returning bogus values, for example)
# >>
# >> I haven't looked at the DCD frame you included yet, but it's important to
# >> know which version of NAMD generated it, since some versions of NAMD
# >> created bad PBC records in the DCD files.  I assume you're using
# >> NAMD 2.5 or later right?  Is it possible you've specified a
# >> set of PBC unit cell basis vectors to NAMD that are not actually
# >> representable with a/b/c/alpha/beta/gamma, but that NAMD is emitting
# >> a gibberish PBC record in the DCD file rather than emitting an error/warning?
# >> Can you send your NAMD config file as well?
# >>
# >> Thanks,
# >>   John
# >>
# >> On Mon, Dec 20, 2004 at 11:36:57AM +0100, Jan Saam wrote:
# >
# >>> > Dear John and other VMD developers,
# >>> >
# >>> > I think I detected a bug in VMD regarding the display of periodic cell
# >>> > images. The error occurs when you load a dcd trajectory that was
# >>> > generated by NAMD using the "dcdUnitCell yes" option.
# >>> > I detected the problem in version 1.8.3a31 while it still worked
# >>> > correctly in VMD 1.8.2
# >>> >
# >>> > It seems like the angles alpha and gamma are assigned in the wrong way.
# >>> > Looking at the source code in Timestep.C:
# >>> > [void Timestep::get_transforms(Matrix4 &a, Matrix4 &b, Matrix4 &c)]
# >>> > shows that the angles alpha and gamma are now defined in a different way
# >>> > as before.
# >>> >
# >>> > before: alpha = angle between A and B
# >>> >         gamma = angle between B and C
# >>> >
# >>> > now:    alpha = angle between B and C
# >>> >         gamma = angle between A and B
# >>> >
# >>> > See the attached files as example:
# >>> > First load the psf/pdb and set the unitcell manually:
# >>> > mol load psf daao_prot_cryw_dows_solv.psf pdb \ daao_prot_cryw_dows_solv.psf
# >>> > # the unitcell according to daao_protfix.xsc:
# >>> > set cell {{60 0 0} {32.1 88 0} {0 0 90}}
# >>> > set_unitcell $cell
# >>> >
# >>> > proc set_unitcell { cell {molid top} } {
# >>> >    # In molinfo the length of the cell vectors and the angles between
# >>> >    # them are saved. $a is assumed to point in x-direction.
# >>> >    set a [veclength [lindex $cell 0]]
# >>> >    set b [veclength [lindex $cell 1]]
# >>> >    set c [veclength [lindex $cell 2]]
# >>> >    set alpha [vecangle [lindex $cell 1] [lindex $cell 2]]
# >>> >    set beta  [vecangle [lindex $cell 0] [lindex $cell 2]]
# >>> >    set gamma [vecangle [lindex $cell 0] [lindex $cell 1]]
# >>> >
# >>> >    molinfo $molid set a $a
# >>> >    molinfo $molid set b $b
# >>> >    molinfo $molid set c $c
# >>> >    molinfo $molid set alpha $alpha
# >>> >    molinfo $molid set beta $beta
# >>> >    molinfo $molid set gamma $gamma
# >>> > }
# >>> >
# >>> > This renders the periodic cell correctly. But now load the dcd file
# >>> > which contains the unitcell info:
# >>> > mol addfile daao_protfix.dcd
# >>> >
# >>> > The dcd contains one frame and the periodic cell is not rendered correctly.
# >>> >
# >>> >
# >>> > Another problem is the way that the periodic is defined in VMD: The
# >>> > vector A is defined to be parallel to the x-axis, while in NAMD I can
# >>> > define a cell where A is arbitrary (e.g. {30 40 0}).
# >>> > These cases will not render correctly in VMD.
# >>> >
# >>> > Merry Christmas,
# >>> >         Jan 

