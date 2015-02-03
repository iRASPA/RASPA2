############################################################
#
# PBCBOX:
# -------
#
# VERSION: 1.2
#
# REQUIREMENTS:
#   The package PBCtools (v1.1 or later) from the VMD script library.
#
# DESCRIPTION:
#     Draws a box around the periodic unit cell boundaries. Also provides a 
#   procedure that automatically updates the box when the frame is changed and 
#   the unit cell changes.  The procedures require the vmd unit cell properties 
#   to be set.
#
# PROCEDURES:
#   pbcbox [OPTIONS...]
#
#     Creates a box that shows the boundaries of the unit cell of 
#   molecule $molid. Returns the list of GIDs.
#
#   Options:
#     -molid $molid     Add the box to the molecule $molid. (default: top)
#     -style lines|dashed|arrows|tubes
#                       Choose the style of the box. (default: lines)
#     -width $w	        Use $width as the width of the lines/arrows/tubes.
#		        (default: 3)
#     -resolution $res  Use $res as the resolution of the cylinders in "tubes" 
#                       style. (default: 8)
#     -origin $origin   $origin has to be a Tcl-list containing three numerical 
#                       values $a $b $c. Place the origin of the box at 
#                       ($a*A, $b*B, $c*C) (default: {-0.5 -0.5 -0.5})
#     -positive         Equal to "-origin {0.0 0.0 0.0}".
#     -parallelepiped|-rectangular
#                       Draw the box be drawn as a parallelpiped,
#                       or as the corresponding rectangular box.
#                       (default: parallelepiped)
#
#
#   pbcbox_update [OPTIONS...]
#
#      Creates a box that shows the boundaries of the unit cell of molecule
#   $molid. The box is automatically redrawn/updated when the frame is changed.
#
#   Options:
#     -state on|off|toggle     
#                       Turn the box on, off or toggle the state.
#     -color $color     Draw the box in $color.
#   All options from the pbcbox procedure can also be used.
#
#  vmd_draw_pbcbox $molid [OPTIONS]
#
#      Procedure to be used with the VMD "draw" procedure. All options from 
#    the pbcbox procedure can be used.
#
# EXAMPLE USAGE:
#   require pbcbox
#   set box [ pbcbox -width 5 -style dashed ]
#   draw delete $box
#
#   pbcbox_update -molid top -style arrows -width 5 -rectangular
#   pbcbox_update
#
#   draw pbcbox -width 7
#
# AUTHOR:
#         Olaf Lenz 
#         olenz _at_ fias.uni-frankfurt.de
#
# This script copies a lot of the ideas and code from Jan Saams script 
# pbctools (that is used anyway) and Axel Kohlmeiers script 
# vmd_draw_unitcell.
#
# Feel free to send comments, bugs, etc.
############################################################

package provide pbcbox 1.2

proc pbcbox { args } {
    # Set the defaults
    set molid       "top"
    set style       "lines"
    set width       3
    set resolution  8
    set positive    0
    set rectangular 0
    set origin_rel  "-0.5 -0.5 -0.5"

    # Parse options
    for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	set arg [ lindex $args $argnum ]
	set val [ lindex $args [expr $argnum + 1]]
	switch -- $arg {
	    "-molid"      { set molid $val; incr argnum }
	    "-style"      { set style $val; incr argnum }
	    "-origin"     { set origin_rel $val; incr argnum }
	    "-positive"   { set origin_rel "0.0 0.0 0.0" }
	    "-parallelepiped" { set rectangular 0 }
	    "-rectangular" { set rectangular 1 }
	    "-width"      { set width $val; incr argnum }
	    "-resolution" { set resolution $val }
	    "-orthogonal" { set rectangular 1 }
	    "-sheared"    {
		switch -- $val {
		    "yes" { set rectangular 0 }
		    "no" { set rectangular 1 }
		    default { set rectangular [expr ! $val] }
		}
		incr argnum
	    }
	    default { puts "unknown option: $arg"; return }
	}
    }

    if { $molid=="top" } then { set molid [ molinfo top ] }

    set cell [ get_unitcell $molid ]
    set A   [lindex $cell 0]
    set B   [lindex $cell 1]
    set C   [lindex $cell 2]
    set origin [lindex $cell 3]
    set origin [vecadd $origin [ vecscale [lindex $origin_rel 0] $A ] ]
    set origin [vecadd $origin [ vecscale [lindex $origin_rel 1] $B ] ]
    set origin [vecadd $origin [ vecscale [lindex $origin_rel 2] $C ] ]

    if { $rectangular } then {
	set A [list [lindex $A 0] 0 0]
	set B [list 0 [lindex $B 1] 0 ]
	set C [list 0 0 [lindex $C 2] ]
    }

    # set up cell vertices
    set vert(0) $origin
    set vert(1) [vecadd $origin $A]
    set vert(2) [vecadd $origin $B]
    set vert(3) [vecadd $origin $A $B]
    set vert(4) [vecadd $origin $C]
    set vert(5) [vecadd $origin $A $C]
    set vert(6) [vecadd $origin $B $C]
    set vert(7) [vecadd $origin $A $B $C]

    set gid ""
    switch $style {
	tubes {
            # set size and radius of spheres and cylinders 
            set srad [expr $width * 0.003 * [veclength [vecadd $A $B $C]]]
            set crad [expr 0.99 * $srad]
	    
            # draw spheres into the vertices ...
            for {set i 0} {$i < 8} {incr i} {
                lappend gid [graphics $molid sphere $vert($i) radius $srad resolution $resolution]
            }
            # ... and connect them with cylinders
            foreach {i j} {0 1  0 2  0 4  1 5  2 3  4 6  1 3  2 6  4 5  7 3  7 5  7 6}  {
                lappend gid [graphics $molid cylinder $vert($i) $vert($j) radius $crad resolution $resolution]
            }
        }
	
        lines {
            set width [expr int($width + 0.5)]
            foreach {i j} {0 1  0 2  0 4  1 5  2 3  4 6  1 3  2 6  4 5  7 3  7 5  7 6}  {
                lappend gid [graphics $molid line $vert($i) $vert($j) width $width style solid]
            }
        }
	
        dashed {
            set width [expr int($width + 0.5)]
            foreach {i j} {0 1  0 2  0 4  1 5  2 3  4 6  1 3  2 6  4 5  7 3  7 5  7 6}  {
                lappend gid [graphics $molid line $vert($i) $vert($j) width $width style dashed]
            }
        }
	
	arrows {
	    set rad [expr $width * 0.003 * [veclength [vecadd $A $B $C]]]
	    foreach { i j } {0 1  0 2  0 4} {
		set middle [vecadd $vert($i) [vecscale 0.9 [vecsub $vert($j) $vert($i) ]]] 
		lappend gid \
		    [graphics $molid cylinder $vert($i) $middle \
			 radius $rad resolution $resolution filled yes ] \
		    [graphics $molid cone $middle $vert($j) \
			 radius [expr $rad * 2.0] resolution $resolution ]
	    }
	}
        default { puts "unknown pbcbox style: $style" ; return }
	
    }

    return $gid
}

# VMD interface for pbcbox (usable with "draw pbcbox")
proc vmd_draw_pbcbox { molid args } {
    return [ eval "pbcbox -molid $molid $args" ]
}

# Create a box that is automatically updated when the frame is changed
# options: 
#  -color (default: blue) 
#  -state (default: toggle)
proc pbcbox_update { args } {
    global pbcbox_gids pbcbox_args pbcbox_color vmd_frame

    # Set the defaults
    set molid "top"
    set state "toggle"
    set color "blue"

    # Parse options
    set pass_args ""
    for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	set arg [ lindex $args $argnum ]
	set val [ lindex $args [expr $argnum + 1]]
	switch -- $arg {
	    "-molid"      { set molid $val; incr argnum }
	    "-color"      { set color $val; incr argnum }
	    "-state"      { set state $val; incr argnum }
	    default { lappend pass_args $arg }
	}
    }

    if { $molid == "top" } then { set molid [ molinfo top ] }

    set pbcbox_color($molid) $color
    set pbcbox_args($molid) $pass_args

    if { [ info exists pbcbox_gids($molid) ] && \
	 ( $state == "off" || $state == "toggle" ) } then {
	puts "Turning off pbcbox."
	display update off
	# deactivate tracing
	trace vdelete vmd_frame($molid) w pbcbox_update_callback
	# delete the box
	pbcbox_update_delete $molid
	display update on
    } elseif { ! [ info exists pbcbox_gids($molid) ] && \
	       ( $state == "on" || $state == "toggle" ) } then {
	puts "Turning on pbcbox:"
	puts "\tcolor=$color"
	puts "\targs=$pass_args"
	display update off
	# activate tracing
	trace variable vmd_frame($molid) w pbcbox_update_callback
	pbcbox_update_draw $molid
	display update on
    }
   
}

# draw the periodic box
proc pbcbox_update_draw { molid } {
    global pbcbox_gids pbcbox_args pbcbox_color
    graphics $molid color $pbcbox_color($molid)
    set pbcbox_gids($molid) \
	[ eval "pbcbox -molid $molid $pbcbox_args($molid)" ] 
}

# delete the periodic box
proc pbcbox_update_delete { molid } {
    global pbcbox_gids
    foreach gid $pbcbox_gids($molid) {
	graphics $molid delete $gid
    }
    array unset pbcbox_gids "$molid"
}

# callback function for vmd_frame, used by "pbcbox update on"
proc pbcbox_update_callback { name1 molid op } {
    display update off
    pbcbox_update_delete $molid
    pbcbox_update_draw $molid
    display update on
}

