getcsize_base<-function(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib,dimerized) {
  #generate edge
  edges<-find_contacts_for_residue(AtomIndexLib, simtraj$coord[,,frame], cutoff_matrix, AtomIndexList, AtomTypeLib,c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame])))
  #format it properly
  edges<-as.matrix(cbind(edges,1))
  edges <- rbind(edges, edges[,c(2,1,3)])
  
  # Add extra edges for dimerized case
  if (dimerized) {
    # Create edges between n and n+1 for odd n (both directions)
    dimer_edges <- matrix(0, nrow = length(AtomIndexLib), ncol = 3)
    edge_count <- 1
    for (n in seq(1, length(AtomIndexLib)-1, by=2)) {
      # Add edge n -> n+1
      dimer_edges[edge_count,] <- c(n, n+1, 1)
      # Add edge n+1 -> n
      dimer_edges[edge_count+1,] <- c(n+1, n, 1)
      edge_count <- edge_count + 2
    }
    # Remove unused rows
    dimer_edges <- dimer_edges[1:(edge_count-1),]
    edges <- rbind(edges, dimer_edges)
  }
  
  attr(edges,"n") <- length(AtomIndexLib)
  #do sna analysis
  compdist<-component.dist(edges)
  #return frame, number of cluster, and largest cluster
  c(frame,max(compdist$csize),length(unique(compdist$membership)))
  
}