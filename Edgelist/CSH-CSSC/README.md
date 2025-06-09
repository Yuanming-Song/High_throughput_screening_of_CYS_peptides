# CSH-CSSC Network Analysis

## Network Statistics

- The file `csh_cssc_network_stats.rda` contains network statistics calculated from the last frame of the entire simulation.
- A new analysis will be performed using the 500 ns frame from the CSSC simulation. The edge list for this frame will be saved as `CSSC_CG_edge_500ns.edges` (and corresponding attribute file `CSSC_CG_edge_500ns.edgeatt`).
- This allows comparison between the network at the end of the simulation and at the 500 ns time point.

## Data and Scripts
- Node and TCL files are managed in the `CG` directory.
- See `copy_cg_nodes_tcl.sh` for automated data transfer. 