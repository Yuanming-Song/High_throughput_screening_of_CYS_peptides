## gyrT_from_edge_tetrapeptide.R

This script provides a function to calculate the gyration tensor moments for the largest cluster in a tetrapeptide system, using a precomputed edgelist (from getEdge). It is designed for use with tetrapeptide simulations where cluster size analysis is handled elsewhere.

**Function:** `gyrT_from_edge_tetrapeptide(edgelist, simtraj, frame, AtomIndexLib)`

- **edgelist**: Edge list for the frame (matrix or data.frame with two columns: node1, node2)
- **simtraj**: Trajectory object (must have $coord and $top)
- **frame**: Frame index (integer)
- **AtomIndexLib**: List mapping residue indices to atom indices

**Returns:**
- List with frame number and gyration tensor moments (from getMoments)

**Note:**
- This function does not perform cluster size analysis. It assumes the largest cluster is of interest.
- The function expects all necessary packages and helper functions (getRogT, getMoments) to be loaded in the environment. 