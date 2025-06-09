# Set library path
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

# Required packages
required_packages <- c("parallel", "doParallel")

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

# Define directories
main_dir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide"
subdata_dir <- file.path(main_dir, "SubData")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript combine_network_stats_tetrapeptide.R <pos> <fromscratch> <single_node> [positions_to_cover]")
}

pos <- as.numeric(args[1])
fromscratch <- as.logical(args[2])
single_node <- as.logical(args[3])
positions_to_cover <- if (length(args) > 3) unique(as.numeric(strsplit(args[4], ",")[[1]])) else 1:4

if (pos < 1 || pos > 4) {
    stop("Error: pos must be an integer between 1 and 4.")
}

# Define residues
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Define states
states <- c("monomer", "dimer")

# Process each position separately
for (curr_pos in positions_to_cover) {
    # combined_results: stores all results for the current position (curr_pos)
    combined_results <- list()
    output_file <- file.path(main_dir, paste0("tetrapeptide_network_stats_C", curr_pos, if(single_node) "_single_node" else "", ".rda"))
    covered <- FALSE
    # If not starting from scratch and output file exists, load previous combined_results from .rda file
    if (!fromscratch && file.exists(output_file)) {
        load(output_file)  # loads 'final_results' from previous run
        # Copy loaded results to combined_results for further merging
        combined_results <- final_results
        covered <- TRUE
    }
    for (state in states) {
        # state_results: accumulates results for this state (monomer/dimer) for the current position
        state_results <- c()
        for (inpres in residues) {
            # Determine expected seqname pattern for this curr_pos and inpres
            if (curr_pos == 1) {
                # CA??
                expected_seqnames <- paste0("C", inpres, rep(residues, each=length(residues)), rep(residues, times=length(residues)))
                expected_seqnames <- paste0("C", inpres, rep(residues, each=length(residues)), rep(residues, times=length(residues)))
                expected_seqnames <- paste0("C", inpres, residues, residues)
                expected_seqnames <- as.vector(outer(residues, residues, function(x, y) paste0("C", inpres, x, y)))
            } else if (curr_pos == 2) {
                # ?CA?
                expected_seqnames <- as.vector(outer(residues, residues, function(x, y) paste0(x, "C", inpres, y)))
            } else if (curr_pos == 3) {
                # ??CA
                expected_seqnames <- as.vector(outer(residues, residues, function(x, y) paste0(x, y, "C", inpres)))
            } else if (curr_pos == 4) {
                # A??C
                expected_seqnames <- as.vector(outer(residues, residues, function(x, y) paste0(inpres, x, y, "C")))
            }
            # Check if all expected seqnames are already present in combined_results[[state]]
            already_covered <- FALSE
            if (!fromscratch && covered && !is.null(combined_results[[state]])) {
                if ("seqname" %in% colnames(combined_results[[state]])) {
                    present_seqnames <- combined_results[[state]][combined_results[[state]]$state == state & substr(combined_results[[state]]$seqname, curr_pos, curr_pos) == inpres, "seqname"]
                    if (length(unique(present_seqnames)) >= length(expected_seqnames)) {
                        already_covered <- TRUE
                    }
                }
            }
            if (!already_covered) {
                # subdata_file: each file contains 'final_results' for a specific residue and state
                subdata_file <- file.path(subdata_dir, 
                                        paste0("Tetrapeptide_network_stats_C", curr_pos, "_", inpres, "_", state, 
                                              if(state == "dimer" && single_node) "_single_node" else "", ".rda"))
                if (file.exists(subdata_file)) {
                    load(subdata_file)  # loads 'final_results' for this residue/state
                    # 'final_results' here is from the subdata file, not the main combined file
                    if (!is.null(final_results[[state]])) {
                        combined_df <- c()
                        for (i in 1:(length(final_results[[state]]))) {
                            combined_df <- rbind(combined_df, final_results[[state]][[i]])
                        }
                        combined_df <- cbind(combined_df, state)
                        state_results <- rbind(state_results, combined_df)
                    }
                    # Remove loaded 'final_results' from subdata file to avoid confusion
                    rm(final_results)
                }
            } else {
                cat(paste0("Skipping subdata for curr_pos=", curr_pos, ", inpres=", inpres, ", state=", state, " (already covered)\n"))
            }
        }
        # Merge with existing results if not starting from scratch and file exists
        if (!fromscratch && covered && !is.null(combined_results[[state]])) {
            combined_results[[state]] <- rbind(combined_results[[state]], state_results)
        } else {
            combined_results[[state]] <- state_results
        }
        # Clear state_results to free memory
        rm(state_results)
    }
    # Save combined results for this position
    final_results <- combined_results  # final_results here is for saving only
    save(final_results, file = output_file)
    cat(paste0("\nCombining complete for position C", curr_pos, "!\n"))
    # Clear variables to avoid duplication in next loop
    rm(combined_results)
    rm(final_results)
}

cat("\nCombining complete!\n") 