getdishis_SIRAH_CSH_CSSC <- function(dimerized,dishis=TRUE) {
  # Initialize variables
  peptide_length <- 1  # CSH/CSSC are single residues
  numatom <- length(which(simtraj$top$resid %in% c("CSX", "CSH") & simtraj$top$resno == 1))
  
  # Example construction of AtomTypeLib
  first_monomer_indices <- 1:numatom
  elety_subset <- simtraj$top$elety[first_monomer_indices]
  
  # Create AtomTypeLib dynamically
  AtomTypeLib <- list()
  unique_elety <- unique(elety_subset)
  
  for (elety in unique_elety) {
    element <- extract_element(elety)
    if (!is.na(element)) {
      AtomTypeLib[[elety]] <- element
    }
  }
  
  # Initialize AtomIndexListPerRes
  AtomIndexListPerRes <- list(Peptide = list(
    BB = NULL
  ))
  
  # Add Sidechains placeholder
  AtomIndexListPerRes$Peptide[["Sidechain1"]] <- NULL
  
  # Select atoms for the current residue
  res_atoms <- simtraj$top[simtraj$top$resno == 1, ]
  
  # Identify backbone atoms (GN GC GO for SIRAH)
  AtomIndexListPerRes$Peptide$BB <- which(res_atoms$elety %in% c("GN", "GC", "GO"))
  
  # Sidechains: All other atoms in the residue
  sidechain_atoms <- which(!res_atoms$elety %in% c("GN", "GC", "GO"))
  AtomIndexListPerRes$Peptide[["Sidechain1"]] <- sidechain_atoms
  
  # Create AtomIndexList
  AtomIndexList <- list()
  for (name in names(AtomTypeLib)) {
    AtomIndexList[[name]] <- which(simtraj$top$elety == name)
  }
  
  # Initialize AtomIndexLib
  AtomIndexLib <- list()
  
  # For CSH-CSSC, we just need to find all molecules with CSX/CSH resid
  csx_indices <- which(simtraj$top$resid %in% c("CSX", "CSH"))
  molecule_count <- 1
  
  # Group atoms into molecules
  while (length(csx_indices) > 0) {
    start_index <- csx_indices[1]
    end_index <- start_index + numatom - 1
    
    # Ensure we have enough atoms for a complete molecule
    if (end_index <= nrow(simtraj$top)) {
      AtomIndexLib[[molecule_count]] <- seq(start_index, end_index)
      molecule_count <- molecule_count + 1
    }
    
    # Remove the indices we just used
    csx_indices <- csx_indices[!(csx_indices %in% seq(start_index, end_index))]
  }
  
  if (dishis) {
    # Get cluster size histogram for all frames
    foreach(frame=framei:dim(simtraj$coord)[3], .combine="+") %dopar% (
      getdishis_base(frame, AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib, dimerized)
    )
  } else {
    foreach(frame=framei:dim(simtraj$coord)[3], .combine="rbind") %dopar% (
      getcsize_base(frame, AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib, dimerized)
    )
  }
} 