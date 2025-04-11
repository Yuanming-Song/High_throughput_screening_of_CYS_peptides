#Generate pdb and psf file of stapled peptide

/dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none -e pfsgen_GCG.tcl -args C G G AA_template_dimer_C3 > psfgen_GCG.log 2>&1

#go to minimization dir, and run it
cd minimization
#load namd
module load namd/2.14b2/gcc.8.4.0
#run min
namd2 min.namd.tcl > min.GCG.log 2>&1 

#go back 
cd ../

#write template pdb
/dfs9/tw/yuanmis1/vmd/bin/vmd -dispdev none C_G_G.psf minimization/C_G_G_min.coor -e getPDBtemplate.tcl > getPDBtemplate.log  2>&1 
