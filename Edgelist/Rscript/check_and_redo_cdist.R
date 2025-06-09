# Load required libraries and source files
.libPaths("/dfs9/tw/yuanmis1/R_libs/")
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

for (packagename in requiredpackages) {
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}

# Source necessary functions
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/generate_sequences_tetrapeptide.R")

# Set directories and parameters
outdir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide"
csizedir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Csizedist/Tetrapeptide/"
binwidth <- 1
poslist <- c(1)#,2,3,4)
# Calculate expected number of sequences
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
expected_cols <- 19^2+ 1  # Number of possible sequences plus the bin column

# Function to process a single sequence
process_sequence <- function(seq, state, single_node) {
  # Load the saved edgelist file
  edgelist_file <- if (state == "dis") {
    file.path(outdir, "dimer", paste0(seq, if (single_node) "_single_node" else "", ".rda"))
  } else {
    file.path(outdir, "monomer", paste0(seq, ".rda"))
  }
  
  if (!file.exists(edgelist_file)) {
    cat(paste0("Missing edgelist file for sequence: ", seq, "\n"))
    return(NULL)
  }
  
  # Load the edgelist
  load(edgelist_file)
  if (!exists("edgelist")) {
    cat(paste0("Invalid edgelist file for sequence: ", seq, "\n"))
    return(NULL)
  }
  
  sizehis <- getdishis_from_edge(edgelist)
  rm(edgelist)
  return(sizehis/sum(sizehis))
}

# Find all cdist files
dimer_files <- list.files(csizedir, pattern="dimer_cdist_Tetrapeptide.*\\.rda$", full.names=TRUE)
monomer_files <- list.files(csizedir, pattern="monomer_cdist_Tetrapeptide\\.rda$", full.names=TRUE)

# Process dimer files
for (dimer_file in dimer_files) {
  # Extract parameters from filename
  filename <- basename(dimer_file)
  parts <- strsplit(filename, "_")[[1]]
  pos <- as.numeric(parts[1])
  inpres <- parts[2]
  single_node <- grepl("single_node", filename)
  
  cat(paste0("Processing file: ", filename, "\n"))
  cat(paste0("Parameters: pos=", pos, ", inpres=", inpres, ", single_node=", single_node, "\n"))
  
  # Generate sequences for this combination
  sequences <- generate_sequences(pos, inpres)
  
  # Load and check the file
  load(dimer_file)
  if (ncol(dimer_sizehis) != expected_cols & pos %in% poslist) {
    cat("Dimer distribution needs fixing...\n")
    # Initialize with bin values
    new_dimer_sizehis <- seq(binwidth, 300+binwidth, binwidth)
    
    # Copy existing valid columns
    if (ncol(dimer_sizehis) > 1) {
      existing_seqs <- colnames(dimer_sizehis)[-1]
      new_dimer_sizehis <- cbind(new_dimer_sizehis, dimer_sizehis[, -1])
    }
    
    # Process missing sequences
    for (seq in sequences) {
      if (!seq %in% colnames(new_dimer_sizehis)) {
        sizehis <- process_sequence(seq, "dis", single_node)
        if (!is.null(sizehis)) {
          new_dimer_sizehis <- cbind(new_dimer_sizehis, sizehis)
          colnames(new_dimer_sizehis)[ncol(new_dimer_sizehis)] <- seq
        } else {
          cat(paste0("Failed to process dimer sequence: ", seq, "\n"))
        }
      }
    }
    
    # Save updated distribution
    dimer_sizehis <- new_dimer_sizehis
    save(dimer_sizehis, file=dimer_file)
  }
}

# Process monomer files
for (monomer_file in monomer_files) {
  # Extract parameters from filename
  filename <- basename(monomer_file)
  parts <- strsplit(filename, "_")[[1]]
  pos <- as.numeric(parts[1])
  inpres <- parts[2]
  
  cat(paste0("Processing file: ", filename, "\n"))
  cat(paste0("Parameters: pos=", pos, ", inpres=", inpres, "\n"))
  
  # Generate sequences for this combination
  sequences <- generate_sequences(pos, inpres)
  
  # Load and check the file
  load(monomer_file)
  if (ncol(monomer_sizehis) != expected_cols) {
    cat("Monomer distribution needs fixing...\n")
    # Initialize with bin values
    new_monomer_sizehis <- seq(binwidth, 300+binwidth, binwidth)
    
    # Copy existing valid columns
    if (ncol(monomer_sizehis) > 1) {
      existing_seqs <- colnames(monomer_sizehis)[-1]
      new_monomer_sizehis <- cbind(new_monomer_sizehis, monomer_sizehis[, -1])
    }
    
    # Process missing sequences
    for (seq in sequences) {
      if (!seq %in% colnames(new_monomer_sizehis)) {
        sizehis <- process_sequence(seq, "mon", FALSE)  # monomer doesn't use single_node
        if (!is.null(sizehis)) {
          new_monomer_sizehis <- cbind(new_monomer_sizehis, sizehis)
          colnames(new_monomer_sizehis)[ncol(new_monomer_sizehis)] <- seq
        } else {
          cat(paste0("Failed to process monomer sequence: ", seq, "\n"))
        }
      }
    }
    
    # Save updated distribution
    monomer_sizehis <- new_monomer_sizehis
    save(monomer_sizehis, file=monomer_file)
  }
} 