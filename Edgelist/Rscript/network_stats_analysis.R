#each monomer as a node?
single_node<-TRUE
# Set library path
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

# Required packages
required_packages <- c("network", "ergm", "parallel", "doParallel", "foreach")

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
            formula <- net ~ edges + kstar(2) + nsp(1) + nsp(2) + esp(0) + esp(1)

            # Calculate network statistics using summary.formula
            stats <- summary(formula)

            return(stats)
        },
        error = function(e) {
            # cat(sprintf("Error calculating network statistics: %s\n", e$message))
            return(c(0, 0, 0, 0, 0, 0))
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
            # Convert stats to dataframe with frame number
            data.frame(
                frame = i,
                value = as.numeric(stats),
                stat_name = names(stats)
            )
        } else {
            # Return empty row with consistent structure
            data.frame(
                frame = i,
                value = 0,
                stat_name = "fail"
            )
        }
    }))

    seqname <- gsub("_", "", seq)
    all_stats <- cbind(all_stats, seqname)

    # Determine C position (1, 2, or 3)
    residues <- unlist(strsplit(seq, "_"))
    c_pos <- which(residues == "C")

    # Get ratio from SASA_files based on state and C position
    ratio_col <- if (state == "monomer") "ratio_mon_fin" else "ratio_dim_fin"
    ratio <- SASA_files[[c_pos]][which(SASA_files[[c_pos]]$label == seqname), ratio_col]
    if (length(ratio) == 0) {
        ratio <- 100
    }
    # Add ratio and C position to stats
    all_stats <- cbind(all_stats, ratio)

    return(all_stats)
}

# Function to process a single sequence
process_single_sequence <- function(seq, state, main_dir) {
    # Get edgelist based on state
    edgelist <- if (state == "monomer") monomer_edgelist[[seq]] else dimer_edgelist[[seq]]

    if (is.null(edgelist)) {
        # Return empty data frame with consistent structure
        empty_df <- data.frame(
            frame = rep(0, 50),
            value = rep(0, 50),
            stat_name = rep("fail", 50),
            seqname = rep(gsub("_", "", seq), 50),
            ratio = rep(0, 50)
        )
        return(empty_df)
    }

    # Calculate network statistics
    stats <- calculate_network_stats_all_frames(edgelist, seq, state)

    # Check if stats is NULL or empty
    if (is.null(stats) || nrow(stats) == 0) {
        # Return empty data frame with consistent structure
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

# Function to generate peptide sequences
generate_peptides <- function(length) {
    all_peptides <- c()

    for (pos in 1:length) {
        if (pos == 1) {
            for (r2 in residues) {
                for (r3 in residues) {
                    all_peptides <- c(all_peptides, paste("C", r2, r3, sep = "_"))
                }
            }
        } else if (pos == 2) {
            for (r1 in residues) {
                for (r3 in residues) {
                    all_peptides <- c(all_peptides, paste(r1, "C", r3, sep = "_"))
                }
            }
        } else if (pos == 3) {
            for (r1 in residues) {
                for (r2 in residues) {
                    all_peptides <- c(all_peptides, paste(r1, r2, "C", sep = "_"))
                }
            }
        }
    }
    return(all_peptides)
}

# Main processing function
# Load SASA files if requested
SASA_files <- list()
for (pos in 1:3) {
    SASA_files[[pos]] <- read.table(
        paste0("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Tripeptide/SASA_score/SASA_result_with_common_mon_ini_C", pos, ".txt"),
        header = FALSE,
        col.names = c("Amino1", "Amino2", "Amino3", "ratio_mon_fin", "ratio_dim_fin")
    )
    SASA_files[[pos]]$label <- paste0(SASA_files[[pos]]$Amino1, SASA_files[[pos]]$Amino2, SASA_files[[pos]]$Amino3)
}

# Define peptide types and their directories
peptide_types <- list(
    tripeptide = list(
        main_dir = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/",
        length = 3
    )
)

# Define states and residues
states <- c("monomer", "dimer")
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Load overall edgelist files
monomer_edgelist <- NULL
dimer_edgelist <- NULL

# Process each peptide type
for (pep_type in names(peptide_types)) {
    cat(sprintf("\nProcessing %s sequences...\n", pep_type))

    # Get directory for this peptide type
    main_dir <- peptide_types[[pep_type]]$main_dir
    pep_length <- peptide_types[[pep_type]]$length

    # Generate sequences for this peptide type
    all_peptides <- generate_peptides(pep_length)
    total_peptides <- length(all_peptides)
    chunk_size <- 500 # Process x peptides at a time
    num_chunks <- ceiling(total_peptides / chunk_size)

    # Initialize results list
    all_results <- list()

    # Process each state and sequence in parallel
    for (state in states) {
        rm(dimer_edgelist)
        rm(monomer_edgelist)
        # Load overall edgelist files
        load(file.path(main_dir, paste0("tripeptide_edgelist_", state, if (single_node & state =="dimer") "_single_node" else "",".rda")))
        assign(paste0(state, "_edgelist"), edgelist)
        cat(sprintf("\nProcessing %s state for %s...\n", state, pep_type))
        rm(edgelist)

        # Process sequences in chunks
        state_results <- list()
        for (chunk_idx in 1:num_chunks) {
            start_idx <- (chunk_idx - 1) * chunk_size + 1
            end_idx <- min(chunk_idx * chunk_size, total_peptides)
            current_chunk <- all_peptides[start_idx:end_idx]

            cat(sprintf(
                "Processing chunk %d/%d (peptides %d-%d)...\n",
                chunk_idx, num_chunks, start_idx, end_idx
            ))

            # Process current chunk in parallel
            chunk_results <- foreach(seq = current_chunk) %dopar% (
                process_single_sequence(seq, state, main_dir)
            )

            # Store chunk results
            state_results[[chunk_idx]] <- chunk_results
            cat(sprintf("Completed chunk %d/%d\n", chunk_idx, num_chunks))
        }

        # Combine all chunk results
        all_results[[state]] <- do.call(c, state_results)
        cat(sprintf("Completed processing %s state\n", state))
    }

    final_results <- all_results
   # do.call(c, all_results)
    # Save results
    save(final_results, file = file.path(main_dir, paste0(pep_type, "_network_stats",if (single_node) "_single_node" else "",".rda")))
    cat(sprintf("\nSaved results to %s\n", file.path(main_dir, paste0(pep_type, "_network_stats.rda"))))
}

cat("\nProcessing complete!\n")

# Example usage:
# process_network_stats(readAP = TRUE)
