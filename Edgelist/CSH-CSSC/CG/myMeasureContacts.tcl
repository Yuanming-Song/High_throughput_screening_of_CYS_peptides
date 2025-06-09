proc myMeasureContacts {sel ref cutoff box } {
	set selpos [$sel get {x y z}]
    set refpos [$ref get {x y z}]
    set selind [$sel get index]
	set refind [$ref get index]
	set selres [$sel get resid]
	set refres [$ref get resid]
	set selList {}
	set refList {}
	foreach vec1 $selpos ind1 $selind res1 $selres {
    	foreach vec2 $refpos ind2 $refind res2 $refres {
    		set dist {}
    		if { $res1 != $res2} {
            	set dif [vecsub $vec1 $vec2]
                foreach x $dif  {
                	lappend dist [expr {$x - $box*round($x/$box)}]
                }
				if {[veclength $dist] <= $cutoff} {
					lappend selList $ind1
					lappend refList $ind2
				}                              
			}
    	}
    }
    return [list $selList $refList]

}

