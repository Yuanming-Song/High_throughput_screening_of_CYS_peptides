gyrT_largest <- function(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib) {
  out<-list()
  edge_list<-find_contacts_for_residue(AtomIndexLib, simtraj$coord[,,frame], cutoff_matrix, AtomIndexList, AtomTypeLib,c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame])))
  #format it properly
  edges<-as.matrix(cbind(edge_list,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  #do sna analysis
  compdist<-component.dist(edges)
  cluster<-which(compdist$membership==which(compdist$csize==max(compdist$csize))[1])
  visited <- c()      # residues whose PBC have been corrected
  queue <- c()        # residues to process next
  box<-c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame]))
  # Start from an arbitrary residue in the cluster, e.g., the first one.
  start <- cluster[1]
  queue <- c(queue, start)
  visited <- c(visited, start)
  clustercoor<-c()
  while(length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    
    # Find all edges that connect current to a neighbor
    # (search in both columns since the graph is undirected)
    neighbors <- unique(c(
      edge_list[edge_list[,1] == current, 2],
      edge_list[edge_list[,2] == current, 1]
    ))
    curcoor<-simtraj$coord[AtomIndexLib[[current]],,frame]
    # 'curcoor' and 'nbcoor' are matrices with 3 columns (x, y, z)
    # 'box' is a numeric vector of length 3, e.g., c(box_x, box_y, box_z)
    
    # Use the first atom in the current residue as the reference coordinate
    ref <- curcoor[1, ]
    
    # Correct PBC for each coordinate (x, y, z)
    for (j in 1:3) {
      # For the current residue:
      diff_cur <- curcoor[, j] - ref[j]
      curcoor[, j] <- curcoor[, j] - round(diff_cur / box[j]) * box[j]
    }
    if(    current ==cluster[1]) {
      clustercoor<-rbind(clustercoor,curcoor)
    }
    # Filter: only consider neighbors that are in the cluster and not yet visited.
    for(nb in neighbors) {
      if(nb %in% cluster && !(nb %in% visited)) {
        nbcoor<-simtraj$coord[AtomIndexLib[[nb]],,frame]
        # Apply the periodic boundary condition correction for neighbor relative to current.
        for (j in 1:3) {
          # For the neighbor residue:
          diff_nb <- nbcoor[, j] - ref[j]
          nbcoor[, j] <- nbcoor[, j] - round(diff_nb / box[j]) * box[j]
        }        
        # Mark neighbor as visited and add it to the queue.
        visited <- c(visited, nb)
        queue <- c(queue, nb)
        clustercoor<-rbind(clustercoor,nbcoor)
        
      }
    }
  }
  #do the gyration tensor thingy
  bigS<-getRogT(clustercoor)
  out[[1]]<-c(frame,max(compdist$csize)[1])
out[[2]]<-getMoments(bigS)
  return(out)
}
