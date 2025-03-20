
set mut_res1 [lindex $argv 0]  ;# Example: LYS
set mut_res2 [lindex $argv 1]  ;# Example: GLY
set mut_res3 [lindex $argv 2]
set mut_res4 [lindex $argv 3]
set template_dir /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/[lindex $argv 4]
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
# Initialize index counter
set index 0
set IndexC -1  ;# Default value if C is not found

# Iterate through each argument in argv
foreach arg $argv {
    # Check if the argument is "C"
    if {$arg == "C"} {
        # Add 1 to the index and save it to a variable
        set IndexC [expr {$index + 1}]
        break  ;# Exit the loop once C is found
    }
    # Increment index
    incr index
}

set res1 [dict get $aa_map $mut_res1]
set res2 [dict get $aa_map $mut_res2]
set res3 [dict get $aa_map $mut_res3]
set res4 [dict get $aa_map $mut_res4]

# Load the topology file
topology /dfs9/tw/yuanmis1/mrsec/ff/toppar_c36_jul22/top_all36_prot.rtf

# Read in the PDB files for the two segments
segment A {
    pdb ${template_dir}/AA_tripeptide_dimer_template_A.pdb
    # Mutate amino acid 1 and 2 of each segment to the specified residues
    mutate 1 $res1
    mutate 2 $res2
    mutate 3 $res3
    mutate 4 $res4
}
regenerate angles dihedrals
# Create a bond between the cysteines at position 3 in segments A and B
coordpdb ${template_dir}/AA_tripeptide_dimer_template_A.pdb A
segment B {
    pdb ${template_dir}/AA_tripeptide_dimer_template_B.pdb
    mutate  1 $res1
    mutate  2 $res2
    mutate 3 $res3
    mutate 4 $res4

}




regenerate angles dihedrals


coordpdb ${template_dir}/AA_tripeptide_dimer_template_B.pdb B

patch DISU A:${IndexC} B:${IndexC}
#standard regeneration for psf file
regenerate angles dihedrals
#guess coordinate for missing atoms (mostly just hydrogen)
guesscoord
# Write the output files
writepdb ${mut_res1}_${mut_res2}_${mut_res3}_${mut_res4}.pdb
writepsf ${mut_res1}_${mut_res2}_${mut_res3}_${mut_res4}.psf
exit