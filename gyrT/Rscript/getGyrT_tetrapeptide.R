# Frame matching: This script matches edgelist frames to simtraj frames by assuming the last edgelist frame matches simtraj$end, and edgelist has step of 1. Frame offset is computed automatically.
# Usage: Rscript getGyrT_tetrapeptide.R <pos> <inpres> [single_node]
single_node<-TRUE
# Set library path
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

#log file for checking simulation length etc
requiredpackages<-c(
  "sna"
  ,"ggplot2"
  ,"doParallel"
  ,"dplyr"
  ,"plotly"
  ,"bio3d"
  ,"geometry"
  ,"bigmemory"
  ,"SOMMD"
  ,"Rcpp"
  ,"tidyr"
)
framei<-200
for (packagename in requiredpackages) {
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}
registerDoParallel(cores = detectCores())

source("/dfs9/tw/yuanmis1/Rscript/fromAlfredo/shape.R")
logfile="csize_anal.log"
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_base.R")
source("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/ML_MD_project_base_R/gyrT_from_edge_tetrapeptide.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/generate_sequences_tetrapeptide.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript getGyrT_tetrapeptide.R <pos> <inpres> [single_node]")
}

pos <- as.numeric(args[1])
inpres <- args[2]
single_node <- if (length(args) >= 3) as.logical(args[3]) else TRUE

if (pos < 1 || pos > 4) {
    stop("Error: pos must be an integer between 1 and 4.")
}

# Set main directory
maindir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI22/"
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y","E","D")
sequences <- generate_sequences(pos, inpres)

# Use 'dis' and 'mon' for state naming
states <- c("dis", "mon")
results <- list()

# Define output directory for gyrT results
outdir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/gyrT/data/Tetrapeptide"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Parallel processing using foreach and %dopar% over sequences
res_list <- foreach(seq = sequences, .packages = c("bio3d", "tidyr")) %dopar% {
  seq_results <- list()
  for (state in states) {
    # Determine edgelist file path
    edgelist_file <- if (state == "dis") {
      file.path("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/dimer", paste0(seq, if (single_node) "_single_node" else "", ".rda"))
    } else {
      file.path("/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/monomer", paste0(seq, ".rda"))
    }
    if (!file.exists(edgelist_file)) {
      seq_results[[state]] <- NULL
      msg <- paste0(seq, " (", state, ") edgelist missing\n")
      cat(msg, file = logfile, append = TRUE)
      next
    }
    load(edgelist_file) # loads 'edgelist'
    # Set directory paths for trajectory
    simdir <- file.path(maindir, paste0("Tetrapeptide_", ifelse(state == "mon", "mon", "dis"), "_C", pos), seq)
    gro_file <- Sys.glob(file.path(simdir, "*md*gro"))
    xtc_file <- Sys.glob(file.path(simdir, "*md*xtc"))
    if (length(gro_file) == 0 || length(xtc_file) == 0) {
      seq_results[[state]] <- NULL
      msg <- paste0(seq, " (", state, ") traj missing\n")
      cat(msg, file = logfile, append = TRUE)
      next
    }
    simtraj <- tryCatch({
      read.trj(xtc_file[1], gro_file[1])
    }, error = function(e) {
      NULL
    })
    if (is.null(simtraj)) {
      seq_results[[state]] <- NULL
      msg <- paste0(seq, " (", state, ") traj load fail\n")
      cat(msg, file = logfile, append = TRUE)
      next
    }
    # AtomIndexLib construction (as in getGyrT)
    resno_sequence <- simtraj$top$resno
    seen_resno <- c()
    numatom <- 0
    for (i in seq_along(resno_sequence)) {
      resno <- resno_sequence[i]
      if (resno %in% seen_resno) {
        if (resno == 1 && length(seen_resno) > 1) {
          numatom <- i - 1
          break
        }
      } else {
        seen_resno <- c(seen_resno, resno)
      }
    }
    peptide_length <- length(seen_resno)
    if (state == "dis") {
      peptide_length <- peptide_length / 2
      numatom <- numatom / 2
    }
    # Set step_size according to single_node logic (see getEdge.R)
    if (single_node) {
      step_size <- numatom
    } else {
      step_size <- if (state == "dis") numatom * 2 else numatom
    }
    AtomIndexLib <- list()
    for (molecule_index in 1:500) {
      start_index <- (molecule_index - 1) * step_size + 1
      end_index <- start_index + step_size - 1
      if (end_index > nrow(simtraj$top)) break
      molecule_block <- simtraj$top[start_index:end_index, ]
      if (all(molecule_block$resid == simtraj$top$resid[1:step_size]) &&
          all(molecule_block$elety == simtraj$top$elety[1:step_size])) {
        AtomIndexLib[[molecule_index]] <- seq(start_index, end_index)
      } else {
        break
      }
    }
    # Frame matching: edgelist last frame matches simtraj$end, step=1
    n_edgelist <- length(edgelist)
    traj_end <- simtraj$end
    frame_offset <- traj_end - n_edgelist + 1
    # For each frame, call gyrT_from_edge_tetrapeptide with correct simtraj frame
    mat <- do.call(rbind, lapply(seq_along(edgelist), function(i) {
      gyrT_from_edge_tetrapeptide(edgelist[[i]], simtraj, i + frame_offset - 1, AtomIndexLib)
    }))
    if (is.null(mat) || nrow(mat) == 0) {
      seq_results[[state]] <- NULL
      next
    }
    colnames(mat) <- c("step", "csize", "T1", "T2", "T3")
    df <- as.data.frame(mat)
    df_long <- pivot_longer(df, cols = c("T1", "T2", "T3"), names_to = "Series", values_to = "Value")
    seq_results[[state]] <- df_long
  }
  list(seq = seq, seq_results = seq_results)
}
gyrT<-list()
# Collect results into gyrT list
for (res in res_list) {
  gyrT[[res$seq]] <- res$seq_results
}

# Save a single gyrT file for all states and sequences for this pos/inpres
outfile <- file.path(outdir, paste0("gyrT_tetrapeptide_", pos, "_", inpres, if (single_node) "_single_node" else "", ".rda"))
save(gyrT, file = outfile)

cat("\nProcessing complete!\n") 