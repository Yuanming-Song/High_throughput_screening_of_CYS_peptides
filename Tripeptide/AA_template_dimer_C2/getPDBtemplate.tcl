#select each segment and write pdb file
foreach seg {A B} {
    set sel [atomselect top "segname $seg"]
    $sel writepdb AA_tripeptide_dimer_template_$seg.pdb
}

exit