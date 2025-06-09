#########################################################################
# TCL script for creating monomer peptide structures using VMD/psfgen.
#
# Usage examples:
#   For tripeptide:
#     vmd -dispdev none -e script.tcl -args A C D
#     -> Creates: A_C_D.pdb and A_C_D.psf
#
#   For pentapeptide:
#     vmd -dispdev none -e script.tcl -args A D C E F
#     -> Creates: A_D_C_E_F.pdb and A_D_C_E_F.psf
#
# Input format:
#   - Takes 2-10 single-letter amino acid codes as arguments
#   - Output files will be named by joining residues with underscore
#   - Example: A B C D -> A_B_C_D.pdb and A_B_C_D.psf
#
# Template selection:
#   - Uses template based on peptide length
#   - Example paths:
#     Tripeptide: .../Tripeptide/AA_template_mon/
#     Pentapeptide: .../Pentapeptide/AA_template_mon/
#
# Note: This is for MONOMER (unstapled) peptides. For stapled peptides,
#       use mutate_template.tcl instead.
#########################################################################

# Get the number of residues
set num_residues [expr {[llength $argv]}]

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

# Construct template directory path for monomer
set template_dir "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/${peptide_type}/AA_template_mon"

# Print information about the peptide being processed
puts "Processing ${peptide_type} monomer"
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

# Process single segment
# Example: for "A C D", will create:
# segment A { pdb template.pdb; mutate 1 ALA; mutate 2 CYS; mutate 3 ASP }
segment A {
    pdb ${template_dir}/AA_template_mon.pdb
    # Mutate each position to the specified residue
    for {set i 0} {$i < $num_residues} {incr i} {
        mutate [expr {$i + 1}] [lindex $residues $i]
    }
}

# Standard regeneration for psf file
regenerate angles dihedrals
coordpdb ${template_dir}/AA_template_mon.pdb A

# Standard regeneration for psf file
regenerate angles dihedrals

# Guess coordinate for missing atoms (mostly just hydrogen)
guesscoord

# Write the output files
# Example: for input A C D, creates A_C_D.pdb and A_C_D.psf
writepdb ${filename}.pdb
writepsf ${filename}.psf
exit