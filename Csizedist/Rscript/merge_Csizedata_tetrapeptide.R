# Define directories and states (maindir is assumed to be defined already)
maindir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Csizedist/"
subdir <- file.path(maindir, "Tetrapeptide/")
states <- c("dimer", "monomer")

# Loop over each state
for (state in states) {
  
  # Create a pattern to match files for the current state; e.g., "^dimer.*\\.rda$"
  pattern <- paste0("*", state, ".*\\.rda$")
  files <- list.files(subdir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found for state:", state, "\n")
    next
  }
  
  # Initialize cumulative data for this state
  state_sizehis_all <- NULL
  
  # Loop over each file for the current state
  for (i in seq_along(files)) {
    load(files[i])
    object_name <- paste0(state, "_sizehis")
    
    if (!exists(object_name)) {
      cat("Object", object_name, "not found in", files[i], "\n")
      next
    }
    
    # For the first file, assign directly; for subsequent files, cbind all columns except the first
    if (i == 1) {
      state_sizehis_all <- get(object_name)
    } else {
      state_sizehis_all <- cbind(state_sizehis_all, get(object_name)[, -1])
    }
  }
  
  # Assign the cumulative result to a variable named <state>_sizehis
  assign(paste0(state, "_sizehis"), state_sizehis_all)
  
  # Save the cumulative clustersizedist to an .rda file under maindir
  out_file <- file.path(maindir, paste0(state, "_cdist_tetrapeptide.rda"))
  save(list = paste0(state, "_sizehis"), file = out_file)
  cat("Saved", state, "clustersizedist to", out_file, "\n")
}