# Set library path
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

# Required packages
required_packages <- c("netdiffuseR", "Matrix","igraph")

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

# Define peptide types and their directories
peptide_types <- list(
  tripeptide = list(
    main_dir = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/",
    base_out_dir = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/adjmat/",
    length = 3
  ),
  tetrapeptide = list(
    main_dir = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/",
    base_out_dir = "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/adjmat/",
    length = 4
  )
)

# Define states
states <- c("monomer", "dimer")

# Define residues excluding "C"
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Function to process edgelist and create adjacency matrix
process_edgelist <- function(edgelist) {
  tryCatch({
    all_nodes <- seq(1,attr(edgelist[[1]],"n"),1) 
    
    # Create a single matrix by combining all frames
    adj_mat <- do.call(cbind, lapply(edgelist, function(edges) {
      g <- graph_from_data_frame(edges, directed = FALSE, vertices = all_nodes)
      mat <- as_adjacency_matrix(g, sparse = FALSE)
      # Ensure consistent row/column order
      mat[all_nodes, all_nodes]
    }))
    
    return(adj_mat)
  }, error = function(e) {
    cat(sprintf("Error processing edgelist: %s\n", e$message))
    return(NULL)
  })
}

# Function to generate peptide sequences
generate_peptides <- function(length) {
  all_peptides <- c()
  
  for (pos in 1:length) {
    if (pos == 1) {
      for (r2 in residues) {
        if (length == 3) {
          for (r3 in residues) {
            all_peptides <- c(all_peptides, paste("C", r2, r3, sep = "_"))
          }
        } else {
          for (r3 in residues) {
            for (r4 in residues) {
              all_peptides <- c(all_peptides, paste("C", r2, r3, r4, sep = "_"))
            }
          }
        }
      }
    } else if (pos == 2) {
      for (r1 in residues) {
        if (length == 3) {
          for (r3 in residues) {
            all_peptides <- c(all_peptides, paste(r1, "C", r3, sep = "_"))
          }
        } else {
          for (r3 in residues) {
            for (r4 in residues) {
              all_peptides <- c(all_peptides, paste(r1, "C", r3, r4, sep = "_"))
            }
          }
        }
      }
    } else if (pos == 3) {
      if (length == 4) {
        for (r1 in residues) {
          for (r2 in residues) {
            for (r4 in residues) {
              all_peptides <- c(all_peptides, paste(r1, r2, "C", r4, sep = "_"))
            }
          }
        }
      }
    } else if (pos == 4) {
      for (r1 in residues) {
        for (r2 in residues) {
          for (r3 in residues) {
            all_peptides <- c(all_peptides, paste(r1, r2, r3, "C", sep = "_"))
          }
        }
      }
    }
  }
  return(all_peptides)
}

# Process each peptide type
for (pep_type in names(peptide_types)) {
  cat(sprintf("\nProcessing %s sequences...\n", pep_type))
  
  # Get directories for this peptide type
  main_dir <- peptide_types[[pep_type]]$main_dir
  base_out_dir <- peptide_types[[pep_type]]$base_out_dir
  pep_length <- peptide_types[[pep_type]]$length
  
  # Create output directories for each state
  for (state in states) {
    dir.create(file.path(base_out_dir, state), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Generate sequences for this peptide type
  all_peptides <- generate_peptides(pep_length)
  
  # Process each state and sequence
  for (state in states) {
    cat(sprintf("\nProcessing %s state for %s...\n", state, pep_type))
    
    # Set state-specific output directory
    state_out_dir <- file.path(base_out_dir, state)
    
    # Process each sequence for this state
    for (seq in all_peptides) {
      #cat(sprintf("Processing sequence %s...\n", seq))
      
      # Load edgelist for this sequence and state
      rda_file <- file.path(main_dir, state, paste0(seq, ".rda"))
      
      if (!file.exists(rda_file)) {
        cat(sprintf("Missing file for sequence %s in %s state\n", seq, state))
        next
      }
      
      # Load the edgelist
      tryCatch({
        out_file <- file.path(state_out_dir, paste0(seq, "_n", total_nodes, "_adjmat.csv"))
        if (file.exists(out_file)) {
          next
        }
        
        load(rda_file)
        
        # Process the edgelist and create adjacency matrices
        adj_mat <- process_edgelist(edgelist)
        
        if (!is.null(adj_mat)) {
          # Get total number of nodes
          total_nodes <- attr(edgelist[[1]], "n")
          
          # Save the adjacency matrix with node count in filename
          write.csv(adj_mat, out_file, row.names = FALSE)
          
        } else {
          cat(sprintf("Failed to create adjacency matrix for sequence %s in %s state\n", seq, state))
        }
      }, error = function(e) {
        cat(sprintf("Error processing sequence %s in %s state: %s\n", seq, state, e$message))
      })
    }
  }
}

cat("\nProcessing complete!\n")

#save it as a matrix 150*(150*202), include frame&num of nodes in file name

#median component size of sequence 

#sufficient statistics 

#sufficient statistics of all sequences as function of time 

