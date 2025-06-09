#
# getDistanceAverage.tcl
#
# JAF
#
# jfreites@uci.edu
#
# version 2 11/20
#
# computes PBCed distances between the COMs of selection pairs
# 
# the output is the mean distance per configuration
#
# This is a VMD TCL script
#


# Input Parameters

# System/Trajectory parameters
#------------------------------------------------------------------------
# dataDir: this is concatenated to the myPSF and file names in theFiles and myReference
# workDir: this is concatenated to the output files
# myPSF: your topology file
# trajFileType: file type for your trajectory file
# step: step length used when reading the trajectory files


set myPSF CSH_CSSC_1to1_50mM_cg_big2_sol_ion.psf
set dataDir ../../csh_cssc_1to1_50mM_cg_big2/
set workDir ""
set step 1
set trajFileType xtc

#-----------------------------------------------------------------------------

# theFiles:
# Provide a TCL list of trajectory file names or use your TCL skills to build it

set theFiles run/CSH_CSSC_1to1_50mM_cg_big2_md.xtc

foreach i {2 3 4 5 6} {
	lappend theFiles run/CSH_CSSC_1to1_50mM_cg_big2_md.part000${i}.xtc
}

# theFileRange:
# Provide a TCL list with the first and last frame number to be analyzed in each
# trajectory file.
# Leave theFileRange empty (set to "") if you want all the frames of all the files
# to be analyzed.

#set theFileRange [list first1 last1 first2 last2 ...]
#

set theFileRange ""

#
#------------------------------------------------------------------------
# Output file name prefix and suffix'
# Full name will include the selection sentences

set outfilepre distances-csh_cssc_1to1_50mM_cg_big2
set oufilesufix part1to6.dat

#set mySelections [list {name S1} {name S1} {name S} ]
#set myReferences [list {name S1} {name S} {name S} ]
#set mySelections [list {name S}  ]
#set myReferences [list {name S}  ]
#set mySelections [list {name BSG}  ]
#set myReferences [list {name BSG}  ]
set mySelections [list {name BSG and resname CSX} {name BSG and resname CSH} {name BSG and resname CSX}  ]
set myReferences [list {name BSG and resname CSX}  {name BSG and resname CSH} {name BSG and resname CSH}]
set cutoff 50.0
#set cutoff 61.0

#-------------------------- END OF User Interface


# Process beg/end trajectory files lists

if {$theFileRange == ""} {
        foreach dcdfile $theFiles {
                append theFileRange "0 -1 "
        }
} else {
        if {[llength $theFileRange] != [expr 2 * [llength $theFiles]]} {
                puts "the File Range list inconsistent with Files list"
                exit
        }
}

set molid [mol new ${dataDir}$myPSF waitfor all]
animate delete all

foreach sele $mySelections refe $myReferences {
	set idf(${sele},${refe}) [open ${workDir}${outfilepre}-[ join ${sele} ""]-[join ${refe} ""]-${oufilesufix} w]
	set sel(${sele},${refe}) [atomselect $molid $sele]
	set ref(${sele},${refe}) [atomselect $molid $refe]
}

set nframes 0
foreach dcdfile $theFiles {begFrame endFrame} $theFileRange {
        animate delete all $molid
        animate read $trajFileType ${dataDir}$dcdfile beg $begFrame end $endFrame skip $step waitfor all $molid
        set totframes [molinfo $molid get numframes]
        for {set currentFrame 0} {$currentFrame < $totframes} {incr currentFrame} {
                incr nframes
                animate goto $currentFrame
                set mybox [molinfo $molid get {a b c}]
                foreach sele $mySelections refe $myReferences {
                	$sel(${sele},${refe}) frame $currentFrame
                	$sel(${sele},${refe}) update
                	$ref(${sele},${refe}) frame $currentFrame
                	$ref(${sele},${refe}) update
                	set selpos [$sel(${sele},${refe}) get {x y z}]
                	set refpos [$ref(${sele},${refe}) get {x y z}]
					set selind [$sel(${sele},${refe}) get index]
					set refind [$ref(${sele},${refe}) get index]
					set selres [$sel(${sele},${refe}) get fragment]
					set refres [$ref(${sele},${refe}) get fragment]
                	set distances {}
                	foreach vec1 $selpos ind1 $selind res1 $selres {
                        foreach vec2 $refpos ind2 $refind res2 $refres {
                                set dist {}
				                if { $res1 != $res2} {
                                	set dif [vecsub $vec1 $vec2]
                                	foreach x $dif box $mybox {
                                	        lappend dist [expr {$x - $box*round($x/$box)}]
                                	}
					                if {[veclength $dist] <= $cutoff} {
                                	#	lappend distances "$ind1 $ind2 [veclength $dist]"
						                 lappend distances [veclength $dist]
				                     }
				                 }
                        }
                	}
					set ave 0.0
					foreach thing $distances {
						set ave [expr $ave + $thing]
					}
                	puts $idf(${sele},${refe}) [expr $ave/[llength $distances] ]
                }
        }
}
foreach sele $mySelections refe $myReferences {
close $idf(${sele},${refe})
}

exit
