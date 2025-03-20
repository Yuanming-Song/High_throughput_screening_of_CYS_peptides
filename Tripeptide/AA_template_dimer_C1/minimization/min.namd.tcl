#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# Min. Warm and Eq.
# PME, Constant Pressure.

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

set Dir /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AA_template_dimer_C1/
set DirFF /dfs9/tw/yuanmis1/mrsec/ff/toppar_c36_jul21/
set sysfile ${Dir}C_G_G
set outfile C_G_G_min
set temperature    300

structure          $sysfile.psf
coordinates        $sysfile.pdb

outputName         $outfile

# Continuing a job from the restart files
if {0} {
    set inputname      npt01
    binCoordinates     $inputname.coor
    binVelocities      $inputname.vel  ;# remove the "temperature" entry if you use this!
    extendedSystem	   $inputname.xsc
}

#firsttimestep      0


#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          ${DirFF}par_all36m_prot.prm

# NOTE: Do not set the initial velocity temperature if you
# have also specified a .vel restart file!
temperature         $temperature


# Periodic Boundary Conditions
# NOTE: Do not set the periodic cell basis if you have also
# specified an .xsc restart file!
if {1} {
    cellBasisVector1     100     0.0     0.0
    cellBasisVector2     0.0       100   0.0
    cellBasisVector3     0.0        0.0    100
    cellOrigin          0.0 0.0 0.0
}


# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.0
switching           on
vdwForceSwitching   yes
switchdist          10.0
pairlistdist        16.0


# Integrator Parameters
timestep            2.0  ;# 1fs/step
rigidBonds          all  ;#
nonbondedFreq       1
fullElectFrequency  2
stepspercycle       20


#PME (for full-system periodic electrostatics)
# values of 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 110, 81, 90, 96,  100, 120, 128
if {1} {
    PME                 yes
    PMEInterpOrder       6
    PMEGridSpacing          1.0
}


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     1     ;# damping coefficient (gamma) of 1/ps
langevinTemp        $temperature

# Constant Pressure Control (variable volume)
if {1} {
    useGroupPressure      yes ;#
    useFlexibleCell       no  ;# no for water box, yes for membrane

    langevinPiston        on
    langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
    langevinPistonPeriod  200.0
    langevinPistonDecay   50.0
    langevinPistonTemp    $temperature
}


restartfreq        10000     ;# 1000steps = every 1ps
restartname ${outfile}
dcdfreq            10000
#xstFreq            10000
outputEnergies     10000
outputTiming	   10000
#outputPressure
#binaryoutput        no


#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################

# Put here any custom parameters that are specific to
# this job (e.g., SMD, TclForces, etc...)
if {0} {
    constraints on
    consexp 2
    consref $sysfile.cons.pdb
    conskfile $sysfile.cons.pdb
    conskcol B
}

#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Minimization and heating
#
minimize 1000
