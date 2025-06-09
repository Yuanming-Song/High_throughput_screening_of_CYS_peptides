#!/usr/bin/env Rscript
.libPaths("/dfs9/tw/yuanmis1/R_libs/")
library(doParallel)
library(foreach)

node_mapping <- function(node) {
  return(floor((node - 1) / 2) + 1)
}

# Function to convert single node edge list to multi node
convert_edge_list <- function(edge_list) {
  # Create new edge list
  new_edges <- list()
  
  # Process each frame's edge list
  for (frame_idx in seq_along(edge_list)) {
    edges <- edge_list[[frame_idx]]
    
    # Convert nodes in edges using the mapping
    # edges[,1] is 'from', edges[,2] is 'to', edges[,3] is weight
    new_edges[[frame_idx]] <- matrix(
      c(
        sapply(edges[,1], node_mapping),
        sapply(edges[,2], node_mapping),
        edges[,3]
      ),
      ncol = 3
    )
    
    # Remove self-edges and duplicates after mapping
    new_edges[[frame_idx]] <- new_edges[[frame_idx]][
      new_edges[[frame_idx]][,1] != new_edges[[frame_idx]][,2],
    ]
    new_edges[[frame_idx]] <- unique(new_edges[[frame_idx]])
    
    # Preserve the 'n' attribute - divide by 2 since we're combining nodes
    attr(new_edges[[frame_idx]], "n") <- attr(edges, "n") / 2
  }
  
  return(new_edges)
}

# Function to process a single file
process_single_file <- function(file) {
  # Load the single_node edge list
  load(file)
  
  # Convert to multi-node version
  new_edgelist <- convert_edge_list(edgelist)
  
  # Create output filename (remove _single_node suffix)
  out_file <- sub("_single_node.rda$", ".rda", file)
  
  # Return both the new edge list and output file path
  return(list(
    edgelist = new_edgelist,
    out_file = out_file
  ))
}

# Function to process files in batches
process_batch <- function(files, batch_size = 19^3) {
  total_files <- length(files)
  num_batches <- ceiling(total_files / batch_size)
  
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, total_files)
    batch_files <- files[start_idx:end_idx]
    
    cat(sprintf("Processing batch %d/%d (files %d-%d)...\n", 
                batch, num_batches, start_idx, end_idx))
    
    # Process files in parallel
    results <- foreach(
      file = batch_files,
      .packages = c("doParallel"),
      .export = c("process_single_file", "convert_edge_list", "node_mapping")
    ) %dopar% {
      process_single_file(file)
    }
    
    # Save results sequentially
    for (result in results) {
      edgelist <- result$edgelist
      save(edgelist, file = result$out_file)
    }
    
    cat(sprintf("Completed batch %d/%d\n", batch, num_batches))
  }
}

# Main function to process all peptide types
process_peptide_type <- function(peptide_type, main_dir) {
  # Define the states to process
  states <- c("dimer", "monomer")
  
  for (state in states) {
    # Skip monomer state as it doesn't have single_node files
    if (state == "monomer") next
    
    # Construct the pattern for single_node files
    pattern <- paste0("*_single_node.rda")
    
    # Get list of all single_node files
    files <- list.files(
      path = file.path(main_dir, peptide_type, state),
      pattern = pattern,
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      cat(sprintf("No single_node files found for %s/%s\n", peptide_type, state))
      next
    }
    
    cat(sprintf("Found %d files for %s/%s\n", length(files), peptide_type, state))
    
    # Process files in batches
    process_batch(files)
  }
}

# Main execution
main_dir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist"

# Set up parallel processing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(sprintf("Using %d cores for parallel processing\n", num_cores))

# Process each peptide type
peptide_types <- c("Dipeptide", "Tripeptide", "Tetrapeptide")
for (type in peptide_types) {
  cat(sprintf("\nProcessing %s...\n", type))
  process_peptide_type(type, main_dir)
}

# Clean up parallel processing
stopCluster(cl) 