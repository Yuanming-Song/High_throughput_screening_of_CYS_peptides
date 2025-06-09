#' Calculate gyration tensor moments for the largest cluster from a given edgelist (tetrapeptide)
#' 
#' @param edgelist Edge list for the frame (matrix or data.frame with two columns: node1, node2)
#' @param simtraj Trajectory object (must have $coord and $top)
#' @param frame Frame index (integer)
#' @param AtomIndexLib List mapping residue indices to atom indices
#' @return Numeric vector: c(frame, max cluster size, T1, T2, T3)
gyrT_from_edge_tetrapeptide <- function(edgelist, simtraj, frame, AtomIndexLib) {
  # Format edge list for sna/component.dist
  edges <- rbind(edgelist, edgelist[, c(2, 1, 3)])
  attr(edges, "n") <- length(AtomIndexLib)
  # Find largest cluster
  compdist <- component.dist(edges)
  max_csize <- max(compdist$csize)
  cluster <- which(compdist$membership == which(compdist$csize == max_csize)[1])
  visited <- c()      # residues whose PBC have been corrected
  queue <- c()        # residues to process next
  # Box dimensions for PBC correction
  box <- c(max(simtraj$coord[, 1, frame]),
           max(simtraj$coord[, 2, frame]),
           max(simtraj$coord[, 3, frame]))
  # Start from an arbitrary residue in the cluster
  start <- cluster[1]
  queue <- c(queue, start)
  visited <- c(visited, start)
  clustercoor <- c()
  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    # Find all edges that connect current to a neighbor
    neighbors <- unique(c(
      edgelist[edgelist[, 1] == current, 2],
      edgelist[edgelist[, 2] == current, 1]
    ))
    curcoor <- simtraj$coord[AtomIndexLib[[current]], , frame]
    ref <- curcoor[1, ]
    # Correct PBC for each coordinate (x, y, z)
    for (j in 1:3) {
      diff_cur <- curcoor[, j] - ref[j]
      curcoor[, j] <- curcoor[, j] - round(diff_cur / box[j]) * box[j]
    }
    if (current == cluster[1]) {
      clustercoor <- rbind(clustercoor, curcoor)
    }
    # Only consider neighbors in the cluster and not yet visited
    for (nb in neighbors) {
      if (nb %in% cluster && !(nb %in% visited)) {
        nbcoor <- simtraj$coord[AtomIndexLib[[nb]], , frame]
        for (j in 1:3) {
          diff_nb <- nbcoor[, j] - ref[j]
          nbcoor[, j] <- nbcoor[, j] - round(diff_nb / box[j]) * box[j]
        }
        visited <- c(visited, nb)
        queue <- c(queue, nb)
        clustercoor <- rbind(clustercoor, nbcoor)
      }
    }
  }
  # Calculate gyration tensor and moments
  bigS <- getRogT(clustercoor)
  moments <- getMoments(bigS)
  # Return as a single row: frame, max cluster size, T1, T2, T3
  return(c(frame, max_csize, moments))
} 