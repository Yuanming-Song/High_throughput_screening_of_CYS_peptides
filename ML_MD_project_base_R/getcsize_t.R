getcsize_base<-function(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib) {
  #generate edge
  edges<-find_contacts_for_residue(AtomIndexLib, simtraj$coord[,,frame], cutoff_matrix, AtomIndexList, AtomTypeLib,c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame])))
  #format it properly
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  #do sna analysis
  compdist<-component.dist(edges)
  #return frame, number of cluster, and largest cluster
  c(frame,max(compdist$csize),length(unique(compdist$membership)))
  
}