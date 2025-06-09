# CSH/CSSC Network Statistics Calculation Script
#
# This script processes edgelists from atomistic (AA) and coarse-grained (CG) CSH and CSSC systems.
# It loads the last frame from each .rda file, regroups the edgelist to the molecule level (each CSH/CSSC is a node),
# adds extra CSSC edges, computes network statistics, and outputs a single data frame with all results.
#
# - For AA: n_nodes = 128
# - For CG CSH: n_nodes = 1284
# - For CG CSSC: n_nodes = 1280
# - n_mol is the number of CSH/CSSC molecules (32, 321, 320, etc.)
# - nodes_per_csh = n_nodes / n_mol
# - For CSSC, an extra edge is added between each odd and next even CSH
#
# Output: A single data frame with all results, saved as an .rda file.

.libPaths("/dfs9/tw/yuanmis1/R_libs/")
library(network)
library(ergm)

# Calculate network statistics for a single frame
calculate_network_stats <- function(edges, n_nodes) {
    tryCatch({
        net <- network(edges[, 1:2],
            directed = FALSE, loops = FALSE, matrix.type = "edgelist",
            num.vertices = n_nodes
        )
        formula <- net ~ edges + kstar(2) + nsp(1) + nsp(2) + esp(0) + esp(1)
        stats <- summary(formula)
        return(stats)
    }, error = function(e) {
        return(c(0, 0, 0, 0, 0, 0))
    })
}

# Explicitly define file paths for each system
csh_aa_file <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/AA/csh32_50mM_every10.edgel.stack.rda"
cssc_aa_file <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/AA/cssc16_50mM_every10.edgel.stack.rda"
csh_cg_file <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/CG/csh_50mM_cg_big_every100_100to407.edgel.stack.rda"
cssc_cg_file <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/CG/cssc_50mM_cg_big_every100_1to359.edgel.stack.rda"

# Regroup edgelist to molecule-level network
# Each CSH (or CSSC) is a single node; edges are between CSH indices
regroup_edgelist <- function(edge_df, nodes_per_csh, n_mol, is_cssc = FALSE) {
    # Map each node to its CSH index
    CSH1 <- ceiling(edge_df[,1] / nodes_per_csh)
    CSH2 <- ceiling(edge_df[,2] / nodes_per_csh)
    edge_df <- data.frame(CSH1, CSH2)
    # Remove self-edges
    edge_df <- edge_df[edge_df$CSH1 != edge_df$CSH2, ]
    # Remove duplicate edges (undirected)
    mol_edges <- t(apply(edge_df, 1, function(x) sort(x)))
    mol_edges <- unique(mol_edges)
    mol_edges <- as.data.frame(mol_edges)
    colnames(mol_edges) <- c("CSH1","CSH2")
    # For CSSC: add extra edge between each odd and next even CSH
    if (is_cssc) {
        extra_edges <- data.frame(
            CSH1 = seq(1, n_mol-1, by=2),
            CSH2 = seq(2, n_mol, by=2)
        )
        mol_edges <- unique(rbind(mol_edges, extra_edges))
    }
    # Add dummy third column for compatibility
    mol_edges <- cbind(mol_edges, 1)
    attr(mol_edges,"n") <- n_mol
    return(mol_edges)
}

# Process last frame from a .rda file, regrouping to molecule-level
process_molecule_level <- function(rda_path, seqname, state, n_mol, is_cssc = FALSE, n_nodes_override = NULL) {
    load(rda_path) # loads 'gs'
    last_frame <- gs[[length(gs)]]
    # Explicitly set n_nodes if override is provided
    if (!is.null(n_nodes_override)) {
        attr(last_frame, "n") <- n_nodes_override
    }
    n_nodes <- attr(last_frame, "n")
    nodes_per_csh <- n_nodes / n_mol
    mol_edges <- regroup_edgelist(last_frame, nodes_per_csh, n_mol, is_cssc)
    # Calculate network statistics and normalize by n_mol
    stats <- calculate_network_stats(mol_edges, n_mol)/n_mol
    stats_df <- data.frame(
        frame = length(gs),
        value = as.numeric(stats),
        stat_name = names(stats),
        seqname = seqname,
        n_nodes = n_nodes,
        nodes_per_csh = nodes_per_csh,
        n_mol = n_mol,
        ratio = NA,
        state = state
    )
    return(stats_df)
}

# Atomistic (AA) systems
csh_aa_stats <- process_molecule_level(csh_aa_file, seqname = "CSH_aa", state = "monomer", n_mol = 32, is_cssc = FALSE, n_nodes_override = 128)
cssc_aa_stats <- process_molecule_level(cssc_aa_file, seqname = "CSH_aa", state = "dimer", n_mol = 32, is_cssc = TRUE, n_nodes_override = 128)

# Coarse-grained (CG) systems
csh_cg_stats <- process_molecule_level(csh_cg_file, seqname = "CSH_cg", state = "monomer", n_mol = 321, is_cssc = FALSE, n_nodes_override = 1284)
cssc_cg_stats <- process_molecule_level(cssc_cg_file, seqname = "CSH_cg", state = "dimer", n_mol = 320, is_cssc = TRUE, n_nodes_override = 1280)

# Combine all results into a single data frame
final_results <- rbind(csh_aa_stats, cssc_aa_stats, csh_cg_stats, cssc_cg_stats)

# Save results
destfile <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats.rda"
save(final_results, file = destfile)
cat("Saved all network statistics to:", destfile, "\n") 