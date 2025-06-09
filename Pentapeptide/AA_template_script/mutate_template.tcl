#########################################################################
# TCL script for creating stapled (dimer) peptide structures using VMD/psfgen.
#
# Usage examples:
#   For tripeptide with C at position 1:
#     vmd -dispdev none -e script.tcl -args C A D
#     -> Creates: C_A_D.pdb and C_A_D.psf
#     -> Uses template: .../Tripeptide/AA_template_dimer_C1/
#
#   For pentapeptide with C at position 3:
#     vmd -dispdev none -e script.tcl -args A D C E F
#     -> Creates: A_D_C_E_F.pdb and A_D_C_E_F.psf
#     -> Uses template: .../Pentapeptide/AA_template_dimer_C3/
#
# Input format:
#   - Takes 2-10 single-letter amino acid codes as arguments
#   - Must include exactly one Cysteine (C) for stapling
#   - Output files will be named by joining residues with underscore
#   - Example: A B C D -> A_B_C_D.pdb and A_B_C_D.psf
#
# Template selection:
#   - Uses template based on peptide length and Cysteine position
#   - Template path format: .../${peptide_type}/AA_template_dimer_C${cys_pos}/
#   - Example: for pentapeptide with C at pos 3 -> .../Pentapeptide/AA_template_dimer_C3/
#
# Stapling:
#   - Creates two identical segments (A and B)
#   - Forms disulfide bond between Cysteine residues of both segments
#   - Example: for sequence A C D, creates:
#     Segment A: ALA-CYS-ASP
#     Segment B: ALA-CYS-ASP
#     Disulfide bond: A:CYS-B:CYS
#########################################################################

# Get the number of residues
set num_residues [expr {[llength $argv]}]

# Find first Cysteine position (1-based index)
# Example: A C D -> cys_pos = 2
set cys_pos -1
for {set i 0} {$i < $num_residues} {incr i} {
    if {[lindex $argv $i] == "C"} {
        set cys_pos [expr {$i + 1}]
        break
    }
}

# Verify Cysteine presence
if {$cys_pos == -1} {
    puts "Error: No Cysteine (C) found in the sequence"
    puts "Note: Stapled peptides require exactly one Cysteine for disulfide bond formation"
    exit 1
}

# Determine peptide type based on length
set peptide_type ""
switch $num_residues {
    2 { set peptide_type "Dipeptide" }
    3 { set peptide_type "Tripeptide" }
    4 { set peptide_type "Tetrapeptide" }
    5 { set peptide_type "Pentapeptide" }
    6 { set peptide_type "Hexapeptide" }
    7 { set peptide_type "Heptapeptide" }
    8 { set peptide_type "Octapeptide" }
    9 { set peptide_type "Nonapeptide" }
    10 { set peptide_type "Decapeptide" }
    default {
        if {$num_residues < 2} {
            puts "Error: Peptide length must be at least 2 residues"
        } else {
            puts "Error: Peptide length must be 10 or fewer residues"
        }
        exit 1
    }
}

# Construct template directory path based on peptide type and Cysteine position
# Example: for pentapeptide with C at pos 3 -> .../Pentapeptide/AA_template_dimer_C3/
set template_dir "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/${peptide_type}/AA_template_dimer_C${cys_pos}"

# Print information about the peptide being processed
puts "Processing ${peptide_type} with Cysteine at position ${cys_pos}"
puts "Using template directory: ${template_dir}"

package require psfgen
pdbalias residue HIS HSD

# Convert single-letter codes to three-letter codes
# Example: A -> ALA, C -> CYS, etc.
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

# Find all positions of Cysteine residues
# Example: A C D -> cys_positions = {2}
set cys_positions {}
for {set i 0} {$i < $num_residues} {incr i} {
    if {[lindex $argv $i] == "C"} {
        lappend cys_positions [expr {$i + 1}]
    }
}

# Convert all residues to three-letter codes
# Example: input A C D -> residues {ALA CYS ASP}
set residues {}
for {set i 0} {$i < $num_residues} {incr i} {
    set one_letter [lindex $argv $i]
    lappend residues [dict get $aa_map $one_letter]
}

# Create filename from all residues
# Example: A C D -> A_C_D
set filename ""
for {set i 0} {$i < $num_residues} {incr i} {
    if {$i > 0} {
        append filename "_"
    }
    append filename [lindex $argv $i]
}

# Load the topology file
topology /dfs9/tw/yuanmis1/mrsec/ff/toppar_c36_jul22/top_all36_prot.rtf

# Process segment A (first copy of the peptide)
# Example: for "A C D", creates:
# segment A { pdb template.pdb; mutate 1 ALA; mutate 2 CYS; mutate 3 ASP }
segment A {
    pdb ${template_dir}/AA_tripeptide_dimer_template_A.pdb
    # Mutate each position to the specified residue
    for {set i 0} {$i < $num_residues} {incr i} {
        mutate [expr {$i + 1}] [lindex $residues $i]
    }
}
regenerate angles dihedrals
coordpdb ${template_dir}/AA_tripeptide_dimer_template_A.pdb A

# Process segment B (second copy of the peptide)
# Creates identical sequence to segment A
segment B {
    pdb ${template_dir}/AA_tripeptide_dimer_template_B.pdb
    # Mutate each position to the specified residue
    for {set i 0} {$i < $num_residues} {incr i} {
        mutate [expr {$i + 1}] [lindex $residues $i]
    }
}
regenerate angles dihedrals
coordpdb ${template_dir}/AA_tripeptide_dimer_template_B.pdb B

# Create disulfide bonds between segments A and B
# Example: for A C D, creates bond between A:2 and B:2 (the Cysteine positions)
foreach cys_pos $cys_positions {
    patch DISU A:$cys_pos B:$cys_pos
}

# Standard regeneration for psf file
regenerate angles dihedrals

# Guess coordinate for missing atoms (mostly just hydrogen)
guesscoord

# Write the output files
# Example: for input A C D, creates A_C_D.pdb and A_C_D.psf
writepdb ${filename}.pdb
writepsf ${filename}.psf
exit