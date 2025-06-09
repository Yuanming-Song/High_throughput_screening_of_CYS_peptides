# want to treat each monomer in dimer as single node?
single_node <- TRUE

.libPaths("/dfs9/tw/yuanmis1/R_libs/")
# log file for checking simulation length etc
requiredpackages <- c(
  "sna",
  "ggplot2",
  "doParallel",
  "dplyr",
  "plotly",
  "bio3d",
  "geometry",
  "bigmemory",
  "SOMMD",
  "Rcpp"
)
statelist <- c("dis", "mon")
framei <- 200
# clustersize bin width
binwidth <- 1
outdir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tripeptide"
csizedir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Csizedist/"

# Load required packages
for (packagename in requiredpackages) {
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}

# Register parallel backend
registerDoParallel(cores = detectCores())

# Source required functions
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getedge_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getEdge.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge.R")

# Set main directory
maindir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide/MARTINI22"

# Function to process a single sequence
process_sequence <- function(seq, pos, maindir, single_node) {
  # Initialize output list
  out <- list(
    output = list(),
    error = character()
  )

  for (state in statelist) {
    # Set directory paths
    simdir <- file.path(maindir, paste0("Tripeptide_", state, "_C", pos), seq)

    # Find simulation files
    gro_file <- Sys.glob(file.path(simdir, "*md*gro"))
    xtc_file <- Sys.glob(file.path(simdir, "*md*xtc"))

    # Check if files exist
    if (length(gro_file) == 0 || length(xtc_file) == 0) {
      out$error <- c(out$error, paste0(seq, " (", state, ") Fail - Missing files\n"))
      next
    }

    # Load trajectory
    simtraj <- tryCatch(
      {
        read.trj(xtc_file[1], gro_file[1])
      },
      error = function(e) {
        out$error <- c(out$error, paste0(seq, " (", state, ") Fail - Trajectory loading error\n"))
        return(out)
        next
      }
    )

    # Check for incomplete trajectories
    if (!is.null(simtraj) && dim(simtraj$coord)[3] < 400) {
      out$error <- c(out$error, paste0(seq, " (", state, ") Incomplete trajectory\n"))
      next
    }

    # Set dimerization state
    dimerized <- ifelse(state == "dis", 1, 0)

    # Get edge list
    edgelist <- getEdge(simtraj, dimerized, single_node)

    # Store results
    if (exists("edgelist")) {
      # Calculate and store sizehis
      sizehis <- getdishis_from_edge(edgelist)
      sizehis <- sizehis / sum(sizehis)
      out$output[[state]] <- list(
        edgelist = edgelist,
        sizehis = sizehis
      )
    }

    rm(edgelist)
  }

  # Add sequence information to output
  out[["sequences"]] <- seq

  # Return the complete output
  return(out)
}

# Main processing loop
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
positions <- c(1, 2, 3)
# Initialize result storage
dimer_sizehis <- seq(binwidth, 300 + binwidth, binwidth)
monomer_sizehis <- seq(binwidth, 300 + binwidth, binwidth)

# Initialize list to store all results
all_results <- list()

for (pos in positions) {
  # Generate sequences for current position
  if (pos == 1) {
    sequences <- outer(residues, residues, function(x, y) paste0("C_", x, "_", y))
  } else if (pos == 2) {
    sequences <- outer(residues, residues, function(x, y) paste0(x, "_C_", y))
  } else {
    sequences <- outer(residues, residues, function(x, y) paste0(x, "_", y, "_C"))
  }

  # Process sequences in parallel
  results <- foreach(seq = 1:length(sequences)) %dopar% {
    process_sequence(sequences[seq], pos, maindir, single_node)
  }

  # Store results for current position
  all_results[[paste0("C", pos)]] <- list(
    results = results,
    sequences = sequences
  )
  # Save all raw results
  save(all_results, file = file.path(outdir, "all_raw_results.rda"))
}


# Process results
for (i in seq_along(results)) {
  seq <- results[[i]]$sequences
  result <- results[[i]]

  # Log errors
  if (length(result$error) > 0) {
    cat(result$error, file = file.path(outdir, "error_log.txt"), append = TRUE)
  }

  # Process output
  for (state in names(result$output)) {
    edgelist <- result$output[[state]]$edgelist

    # Save edgelist
    if (state == "dis") {
      save(edgelist, file = file.path(outdir, "dimer", paste0(seq, if (single_node) "_single_node" else "", ".rda")))
    } else {
      save(edgelist, file = file.path(outdir, "monomer", paste0(seq, ".rda")))
    }

    sizehis <- result$output[[state]]$sizehis

    if (state == "dis") {
      dimer_sizehis <- cbind(dimer_sizehis, sizehis)
      colnames(dimer_sizehis)[ncol(dimer_sizehis)] <- seq
    } else {
      monomer_sizehis <- cbind(monomer_sizehis, sizehis)
      colnames(monomer_sizehis)[ncol(monomer_sizehis)] <- seq
    }
  }
}
# Save final sizehis results
save(dimer_sizehis, file = file.path(csizedir, paste0("dimer_cdist_tripeptide", if (single_node) "_single_node" else "", ".rda")))
save(monomer_sizehis, file = file.path(csizedir, "monomer_cdist_tripeptide.rda"))
