

set template_dir /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Pentapeptide/ 
package require psfgen
pdbalias residue HIS HSD
# Specify the mutations using single-letter amino acid codes

# Convert single-letter codes to three-letter codes
set aa_map {
    A ALA
    R ARG
    N ASN
    D ASP
    C CYS
    E GLU
    Q GLN
    G GLY
    H HSD
    I ILE
    L LEU
    K LYS
    M MET
    F PHE
    P PRO
    S SER
    T THR
    W TRP
    Y TYR
    V VAL
}
# Assume the script is called with one argument: the C index (an integer 1–4)
# If the provided index corresponds to a position, that position gets "CYS" and all others "PHE"

# Check that exactly one argument is provided
if {[llength $argv] != 1} {
    puts "Usage: script.tcl <C_index (1-4)>"
    exit 1
}

set c_index [lindex $argv 0]

# Validate that c_index is an integer between 1 and 4
if {![string is integer -strict $c_index] || $c_index < 1 || $c_index > 4} {
    puts "Error: C index must be an integer between 1 and 4."
    exit 1
}

# Default: all residues are "PHE"
set res_list [list "PHE" "PHE" "PHE" "PHE"]

# Set the residue at the specified index (adjust for 0-indexing) to "CYS"
lset res_list [expr {$c_index - 1}] "CYS"

# Extract individual residues
set res1 [lindex $res_list 0]
set res2 [lindex $res_list 1]
set res3 [lindex $res_list 2]
set res4 [lindex $res_list 3]
# Load the topology file
topology /dfs9/tw/yuanmis1/mrsec/ff/toppar_c36_jul22/top_all36_prot.rtf

# Read in the PDB files for the two segments
segment A {
    pdb ${template_dir}/FFFF_origin.pdb
    # Mutate amino acid 1 and 2 of each segment to the specified residues
    mutate 1 $res1
    mutate 2 $res2
    mutate 3 $res3
        mutate 4 $res4

}
regenerate angles dihedrals
# Create a bond between the cysteines at position 3 in segments A and B
coordpdb ${template_dir}/FFFF_origin.pdb A
segment B {
    pdb ${template_dir}/FFFF_moved.pdb
    mutate  1 $res1
    mutate  2 $res2
    mutate 3 $res3
            mutate 4 $res4

}




regenerate angles dihedrals


coordpdb ${template_dir}/FFFF_moved.pdb B

patch DISU A:${c_index} B:${c_index}
#standard regeneration for psf file
regenerate angles dihedrals
#guess coordinate for missing atoms (mostly just hydrogen)
guesscoord
# Write the output files
writepdb Template_dimer.pdb
writepsf Template_dimer.psf
exit