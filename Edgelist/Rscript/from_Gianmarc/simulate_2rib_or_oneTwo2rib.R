# This code demonstrates how to run ergm simulations of 2-ribbons and 1,2 2-ribbons.
# It also uses the fibril_assay function to determine which nodes belong to one
# of these two fibril topologies, and also returns the fibril yeild for each type.
# 
# Gianmarc Grazioli, Ph.D.
# 
# NOTE: You must change the working directory below to your working directory:
setwd("~/Desktop/research/projects/polymorphic_ERGM/starterCodes/v3")

library(ergm)
library(sna)
library(parallel)
source("create_target_structure.R")
source("fibril_assay.R")

begin=proc.time()[[3]]
# set.seed(301656) # set random seed if consistent results desired

nCount=48 #number of nodes

# 2-ribbon tau: 
# tau = c(109 - log(nCount), -24.799511, -1.121684,  2.953991, 0.000000)

# 1,2 2-ribbon tau: 
tau=c(157.912 -log(nCount),-24,-3.3,-2.7,-20)

simCount<-5 # number of times to repeat simulation

# MCMC burn time:
t<-sum(2^(1:20)) # later change to (1:22)

np=2 #number of processors to use

sims<-vector(mode="list",length=simCount)
initGraph<-network.initialize(nCount,directed=F)
for(i in 1:simCount){
  cat("Working on rep=",i,"\n")
  sims[[i]]<-simulate(initGraph~edges+kstar(2)+nsp(1:2)+esp(0), coef=tau, control=control.simulate.formula(MCMC.burn=t,MCMC.interval=1),constraint=~bd(maxout=12),mc.cores=np)
  #gplot3d(sims[[i]],gmode="graph")
  plot(sims[[i]])
}

# See which nodes are part of a fibril structure:
assayList <- list()
for (i in 1:length(sims)) {
  assayList[[i]] <- fibril_assay(sims[[i]])
}

fibrilYields <- list()
for (i in 1:length(assayList)){
  fibrilYields[[i]] <- colSums(assayList[[i]])/nrow(assayList[[i]])
}

fibrilYields

# Added an example of counting sufficent statistics for Song here:
statListExample <- summary(sims[[2]] ~ edges + kstar(2) + nsp(1:2) + esp(0))
# It returns a vector of the chosen sufficient statistics for the graph stored at sims[[2]]

statListExample

end=proc.time()[[3]]
cat("Calculation took ", end-begin, "seconds.\n")
