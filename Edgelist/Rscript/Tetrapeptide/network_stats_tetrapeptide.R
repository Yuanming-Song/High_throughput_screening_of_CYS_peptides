# Set library path
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

# Required packages
required_packages <- c("network", "ergm", "parallel", "doParallel", "foreach")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript network_stats_tetrapeptide.R <pos> <inpres> [single_node]")
}

pos <- as.numeric(args[1])
inpres <- args[2]
single_node <- if (length(args) >= 3) as.logical(args[3]) else TRUE

if (pos < 1 || pos > 4) {
    stop("Error: pos must be an integer between 1 and 4.")
}

# Function to safely load packages
load_packages <- function(packages) {
    for (package in packages) {
        if (!requireNamespace(package, quietly = TRUE)) {
            install.packages(package)
            library(package, character.only = TRUE)
        } else {
            library(package, character.only = TRUE)
        }
    }
}

# Load required packages
load_packages(required_packages)
registerDoParallel(cores = detectCores())

# Function to calculate network statistics for a single frame
calculate_network_stats <- function(edges, n_nodes) {
    tryCatch(
        {
            # Create network object
            net <- network(edges[, 1:2],
                directed = FALSE, loops = FALSE, matrix.type = "edgelist",
                num.vertices = n_nodes
            )

            # Create formula for statistics
            formula <- net ~ edges + kstar(2) + nsp(1) + nsp(2) + esp(0) + esp(1) + cycle(4)

            # Calculate network statistics using summary.formula
            stats <- summary(formula)

            return(stats)
        },
        error = function(e) {
            return(c(0, 0, 0, 0, 0, 0, 0))
        }
    )
}

# Function to calculate network statistics for all frames
calculate_network_stats_all_frames <- function(edgelist, seq, state) {
    n_nodes <- attr(edgelist[[1]], "n")

    # Process all frames in edgelist
    all_stats <- do.call(rbind, lapply(seq_along(edgelist), function(i) {
        stats <- calculate_network_stats(edgelist[[i]], n_nodes)
        if (!is.null(stats)) {
            data.frame(
                frame = i,
                value = as.numeric(stats),
                stat_name = names(stats)
            )
        } else {
            data.frame(
                frame = i,
                value = 0,
                stat_name = "fail"
            )
        }
    }))

    seqname <- gsub("_", "", seq)
    all_stats <- cbind(all_stats, seqname)

    # Determine C position (1, 2, 3, or 4)
    residues <- unlist(strsplit(seq, "_"))
    c_pos <- which(residues == "C")

    # Get ratio from SASA_files based on state and C position
    ratio_col <- if (state == "monomer") "ratio_mon_fin" else "ratio_dim_fin"
    ratio <- SASA_files[[c_pos]][which(SASA_files[[c_pos]]$label == seqname), ratio_col]
    if (length(ratio) == 0) {
        ratio <- 100
    }
    all_stats <- cbind(all_stats, ratio)

    return(all_stats)
}

# Function to process a single sequence
process_single_sequence <- function(seq, state, main_dir) {
    edgelist <- if (state == "monomer") monomer_edgelist[[seq]] else dimer_edgelist[[seq]]

    if (is.null(edgelist)) {
        empty_df <- data.frame(
            frame = rep(0, 50),
            value = rep(0, 50),
            stat_name = rep("fail", 50),
            seqname = rep(gsub("_", "", seq), 50),
            ratio = rep(0, 50)
        )
        return(empty_df)
    }

    stats <- calculate_network_stats_all_frames(edgelist, seq, state)

    if (is.null(stats) || nrow(stats) == 0) {
        empty_df <- data.frame(
            frame = rep(0, 50),
            value = rep(0, 50),
            stat_name = rep("fail", 50),
            seqname = rep(gsub("_", "", seq), 50),
            ratio = rep(0, 50)
        )
        return(empty_df)
    }

    return(stats)
}

# Function to generate tetrapeptide sequences for a specific position and first residue
generate_sequences <- function(pos, inpres) {
    residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
    sequences <- c()
    
    if (pos == 1) {
        for (r2 in residues) {
            for (r3 in residues) {
                for (r4 in residues) {
                    if (r2 == inpres) {
                        sequences <- c(sequences, paste("C", r2, r3, r4, sep = "_"))
                    }
                }
            }
        }
    } else if (pos == 2) {
        for (r1 in residues) {
            for (r3 in residues) {
                for (r4 in residues) {
                    if (r1 == inpres) {
                        sequences <- c(sequences, paste(r1, "C", r3, r4, sep = "_"))
                    }
                }
            }
        }
    } else if (pos == 3) {
        for (r1 in residues) {
            for (r2 in residues) {
                for (r4 in residues) {
                    if (r1 == inpres) {
                        sequences <- c(sequences, paste(r1, r2, "C", r4, sep = "_"))
                    }
                }
            }
        }
    } else if (pos == 4) {
        for (r1 in residues) {
            for (r2 in residues) {
                for (r3 in residues) {
                    if (r1 == inpres) {
                        sequences <- c(sequences, paste(r1, r2, r3, "C", sep = "_"))
                    }
                }
            }
        }
    }
    return(sequences)
}

# Main processing
# Load SASA files
SASA_files <- list()
for (pos_idx in 1:4) {
    SASA_files[[pos_idx]] <- read.table(
        paste0("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Tetrapeptide/SASA_score/SASA_result_with_common_mon_ini_tetra_C", pos_idx, ".txt"),
        header = FALSE,
        col.names = c("Amino1", "Amino2", "Amino3", "Amino4", "ratio_mon_fin", "ratio_dim_fin")
    )
    SASA_files[[pos_idx]]$label <- paste0(SASA_files[[pos_idx]]$Amino1, SASA_files[[pos_idx]]$Amino2, SASA_files[[pos_idx]]$Amino3, SASA_files[[pos_idx]]$Amino4)
}

# Define directories
main_dir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide"
subdata_dir <- file.path(main_dir, "SubData")

# Define states
states <- c("monomer", "dimer")

# Generate sequences for this position and inpres
sequences <- generate_sequences(pos, inpres)

# Process each state
for (state in states) {
    # Load edgelist files
    edgelist_file <- file.path(subdata_dir, paste0("Tetrapeptide_edgelist_", state, "_C", pos, "_", inpres,
                              if(state == "dimer" && single_node) "_single_node" else "", ".rda"))
    load(edgelist_file)
    assign(paste0(state, "_edgelist"), edgelist)
    
    # Check for missing sequences in edgelist
    edgelist_names <- names(edgelist)
    missing_sequences <- setdiff(sequences, edgelist_names)
    if (length(missing_sequences) > 0) {
        log_file <- file.path(main_dir, "network_stats_missing.log")
        log_lines <- paste(missing_sequences, state, sep = "\t")
        # If file doesn't exist, create it with header
        if (!file.exists(log_file)) {
            writeLines("sequence\tstate", log_file)
        }
        write(log_lines, file = log_file, append = TRUE)
    }
    rm(edgelist)
    
    # Process sequences
    results <- foreach(seq = sequences) %dopar% {
        process_single_sequence(seq, state, main_dir)
    }
    
    # Save results for this state and inpres
    final_results <- list()
    final_results[[state]] <- results
    output_file <- file.path(subdata_dir, 
                            paste0("Tetrapeptide_network_stats_C", pos, "_", inpres, "_", state,
                                  if(state == "dimer" && single_node) "_single_node" else "", ".rda"))
    save(final_results, file = output_file)
}

cat("\nProcessing complete!\n") 