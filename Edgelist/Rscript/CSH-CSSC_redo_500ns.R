# CSH/CSSC Network Statistics Calculation Script (500ns CSSC)
#
# This script loads csh_cssc_network_stats.rda for CSH (monomer) info,
# but for CSSC (dimer), reads the newly generated 500ns edges file,
# processes it as a molecule-level network, and outputs a combined data frame.
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

# Regroup edgelist to molecule-level network
regroup_edgelist <- function(edge_df, nodes_per_csh, n_mol, is_cssc = FALSE) {
    CSH1 <- ceiling(edge_df[,1] / nodes_per_csh)
    CSH2 <- ceiling(edge_df[,2] / nodes_per_csh)
    edge_df <- data.frame(CSH1, CSH2)
    edge_df <- edge_df[edge_df$CSH1 != edge_df$CSH2, ]
    mol_edges <- t(apply(edge_df, 1, function(x) sort(x)))
    mol_edges <- unique(mol_edges)
    mol_edges <- as.data.frame(mol_edges)
    colnames(mol_edges) <- c("CSH1","CSH2")
    if (is_cssc) {
        extra_edges <- data.frame(
            CSH1 = seq(1, n_mol-1, by=2),
            CSH2 = seq(2, n_mol, by=2)
        )
        mol_edges <- unique(rbind(mol_edges, extra_edges))
    }
    mol_edges <- cbind(mol_edges, 1)
    attr(mol_edges,"n") <- n_mol
    return(mol_edges)
}

# Load all results from rda
load("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats.rda")

# Extract all but CSH_cg dimer (i.e., keep CSH_aa monomer, CSH_aa dimer, CSH_cg monomer)
keep_rows <- !(final_results$seqname == "CSH_cg" & final_results$state == "dimer")
final_results_base <- final_results[keep_rows, ]

# Read the 500ns CSSC edge file (for CSH_cg dimer, 500ns)
cssc_edges <- read.table("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/CG/CSSC_CG_edge_500ns.edges", header=FALSE)
n_mol_cssc <- 320
n_nodes_cssc <- 1280
nodes_per_csh_cssc <- n_nodes_cssc / n_mol_cssc
cssc_mol_edges <- regroup_edgelist(cssc_edges, nodes_per_csh_cssc, n_mol_cssc, is_cssc=TRUE)
cssc_stats <- calculate_network_stats(cssc_mol_edges, n_mol_cssc)/n_mol_cssc
cssc_stats_df <- data.frame(
    frame = 500,
    value = as.numeric(cssc_stats),
    stat_name = names(cssc_stats),
    seqname = "CSH_cg",
    n_nodes = n_nodes_cssc,
    nodes_per_csh = nodes_per_csh_cssc,
    n_mol = n_mol_cssc,
    ratio = NA,
    state = "dimer"
)

# Combine all results
final_results_500ns <- rbind(final_results_base, cssc_stats_df)

# Save results
save(final_results_500ns, file = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats_500ns.rda")
cat("Saved 500ns network statistics to: /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats_500ns.rda\n") 