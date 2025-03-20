# Get the argument peplength
set peplength 3
#[lindex $argv 0]

# Loop through each residue in chain B from 1 to peplength
for {set resid 1} {$resid <= $peplength} {incr resid} {
    # Select the residue with chain B and increment its resid
    set sel [atomselect top "resid $resid and chain B"]
    $sel set resid [expr {$resid + $peplength}]
    $sel delete
}

# Select all atoms and change their chain to A
set all_sel [atomselect top "all"]
$all_sel  writepdb aa_for_SIRAH.pdb
exit