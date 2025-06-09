#
# getEdgesSmallsirah.tcl (single-GRO version)
#
# This script reads one GRO file (500ns frame) and generates an edgelist for that structure.
# Input: GRO file (full path), node list file (in current dir)
# Output: CSSC_CG_edge_500ns.edges and CSSC_CG_edge_500ns.edgeatt

# Set the node list file
set myNodeList CSSC_50mM_cg_sol_ion.psf.nodes

# Set the full path to the 500ns GRO file
set grofile "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/CSH-CSSC/CSSC_50mM_cg_md_500ns.gro"

# Load the structure
set molwork [mol new $grofile]

# Set the output file name
set outfile CSSC_CG_edge_500ns

# Source the contact measurement script (must be in current dir)
source myMeasureContacts.tcl

# Main selection
set mySelection "resname CSX CSH"

# Cutoff Parameters
set Cutoff 5.0
set carbonCutoff $Cutoff
set sulfurCutoff $Cutoff

# Create a dictionary based on Eric's node list
set idf [open $myNodeList r]
while {[gets $idf line] >= 0} {
    set nindex [lindex $line 0]
    dict set network $nindex resname [lindex $line 1]
    dict set network $nindex resid [lindex $line 2]
    dict set network $nindex moid [lindex $line 3]
    dict set network $nindex nname [lindex $line 4]
    dict set network $nindex ntype [lindex $line 5]
    foreach thing [lrange $line 6 end] {
        dict set indices $thing $nindex
    }
}
close $idf

# Atom selections
set nohs [atomselect $molwork "($mySelection)"]
set nocarbons [atomselect $molwork "($mySelection) and not name BCG BCE1 BCE2 GC2 BCG"]
set carbonAtoms [atomselect $molwork "($mySelection) and name BCG BCE1 BCE2 GC2 BCG"]
set sulfurAtoms [atomselect $molwork "($mySelection) and name BSG BPG"]

# Only one frame, analyze directly
set idf [open ${outfile}.edges w]
set idf2 [open ${outfile}.edgeatt w]

$carbonAtoms update
$nocarbons update
$nohs update
set box [molinfo $molwork get a]
set carbonContacts [myMeasureContacts $carbonAtoms $carbonAtoms $carbonCutoff $box]
set sulfurContacts [myMeasureContacts $sulfurAtoms $sulfurAtoms $sulfurCutoff $box]
set Contacts [myMeasureContacts $nocarbons $nohs $Cutoff $box]
set Ai [concat [lindex $carbonContacts 0] [lindex $Contacts 0] [lindex $sulfurContacts 0]]
set Aj [concat [lindex $carbonContacts 1] [lindex $Contacts 1] [lindex $sulfurContacts 1]]
set weight {}
set interactType {}
foreach thing1 $Ai thing2 $Aj {
    set Ni [dict get $indices $thing1]
    set Nj [dict get $indices $thing2]
    if {$Ni != $Nj} {
        set pair [lsort -integer [list $Ni $Nj]]
        dict incr weight $pair
        if {![dict exists $interactType $pair]} {
            set interaction [list [dict get $network [lindex $pair 0] nname] [dict get $network [lindex $pair 1] nname]]
            switch -- $interaction {
                {PHE PHE} {
                    dict set interactType $pair PHE-PHE
                }
                {AM1 AM2} -
                {AM2 AM1} -
                {AM1 AM1} -
                {AM2 AM2} {
                    dict set interactType $pair AM-AM
                }
                {CYS CYS} {
                    dict set interactType $pair DISU
                }
                {AM1 PHE} -
                {PHE AM1} {
                    dict set interactType $pair AM1-PHE
                }
                {AM2 PHE} -
                {PHE AM2} {
                    dict set interactType $pair AM2-PHE
                }
                default {
                    dict set interactType $pair STER
                }
            }
        }
    }
}
dict for {key val} $weight {
    puts $idf "$key $val"
}
dict for {key val} $interactType {
    puts $idf2 "$key $val"
}
close $idf
close $idf2
puts "edge list written for single structure"

exit
