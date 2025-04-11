getedge_base<-function(simtraj,frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib, single_node=FALSE) {
  #generate edge
  edges<-find_contacts_for_residue(AtomIndexLib, simtraj$coord[,,frame], cutoff_matrix, AtomIndexList, AtomTypeLib,c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame])))
  #format it properly
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  if (single_node) {
    #make sure disulfide bond all counted
    edgetemp<-cbind(seq(1,length(AtomIndexLib)-1,2),seq(2,length(AtomIndexLib),2),1)
    edgetemp<-rbind(edgetemp,edgetemp[,c(2,1,3)])
    
    # Fast duplicate checking using matrix operations
    is_duplicate <- duplicated(rbind(edges[,1:2], edgetemp[,1:2]))
    new_edges <- edgetemp[!is_duplicate[(nrow(edges)+1):length(is_duplicate)],]
    if(nrow(new_edges) > 0) {
      edges <- rbind(edges, new_edges)
    }
  }
    attr(edges,"n")<-length(AtomIndexLib)
  #return edge list
  edges
}
