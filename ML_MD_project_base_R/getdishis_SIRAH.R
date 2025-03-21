getdishis_SIRAH<-function(dimerized,dishis=TRUE) {
  
  # Initialize variables
  resno_sequence <- simtraj$top$resno
  seen_resno <- c()    # Track seen residues
  numatom <- 0         # Number of rows in the first peptide monomer
  
  # Analyze line-by-line until resno starts over after full coverage
  for (i in seq_along(resno_sequence)) {
    resno <- resno_sequence[i]
    
    if (resno %in% seen_resno) {
      # If resno == 1 and all other residues are seen, we stop
      if (resno == 1 && length(seen_resno) > 1) {
        numatom <- i - 1  # Rows covering the first monomer
        break
      }
    } else {
      seen_resno <- c(seen_resno, resno)  # Add unseen residue
    }
  }
  
  # Determine peptide length
  peptide_length <- length(seen_resno)
  if (dimerized) {
    peptide_length<-peptide_length/2
    numatom<-numatom/2
  }
  
  
  
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
  
  last_residue <- peptide_length
  
  # Initialize AtomIndexListPerRes
  AtomIndexListPerRes <- list(Peptide = list(
    BB = NULL
  ))
  
  # Add Sidechains placeholders for each residue
  for (i in 1:peptide_length) {
    AtomIndexListPerRes$Peptide[[paste0("Sidechain", i)]] <- NULL
  }
  
  # Loop through each residue
  for (trures_index in 1:length(seen_resno)) {
    if (trures_index>peptide_length) {
      res_index<-trures_index-peptide_length
    } else {
      res_index<-trures_index 
    }
    
    # Select atoms for the current residue
    res_atoms <- simtraj$top[simtraj$top$resno == trures_index, ]
    
    # Identify Nterm (GN in first residue)
    #need to be GN GC GO for SIRAH
    AtomIndexListPerRes$Peptide$BB <- c(AtomIndexListPerRes$Peptide$BB,which(res_atoms$elety %in% c("GN", "GC", "GO")))
    
    
    # Sidechains: All other atoms in the residue
    sidechain_atoms <- which(!res_atoms$elety %in% c("GN", "GC", "GO") )
    AtomIndexListPerRes$Peptide[[paste0("Sidechain", res_index)]] <- sidechain_atoms
  }
  AtomIndexList<-list()
  for (name in names(AtomTypeLib)) {
    AtomIndexList[[name]]<-which(simtraj$top$elety==name) 
  }
  # Initialize AtomIndexLib
  AtomIndexLib <- list()
  
  
  # Determine step size
  step_size <- if (dimerized) numatom * 2 else numatom
  
  # Iterate through `resno` for molecules
  for (molecule_index in 1:2000) {
    # Calculate start and end indices for this molecule
    start_index <- (molecule_index - 1) * step_size + 1
    end_index <- start_index + step_size - 1
    
    # Ensure we are within bounds of `simtraj$top`
    if (end_index > nrow(simtraj$top)) break
    
    # Check if this slice matches the repeating pattern
    molecule_block <- simtraj$top[start_index:end_index, ]
    if (
      all(molecule_block$resid == simtraj$top$resid[1:step_size]) &&
      all(molecule_block$elety == simtraj$top$elety[1:step_size]) ) {
      # Record line numbers for this molecule
      AtomIndexLib[[molecule_index]] <- seq(start_index, end_index)
    } else {
      break
    }
  }
  if (dishis) {
    #get cluster size histogram for all frame, and add them up, return the result directly
    foreach(frame=framei:dim(simtraj$coord)[3],.combine="+") %dopar% (
      getdishis_base(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib)
    )
  } else {
    foreach(frame=framei:dim(simtraj$coord)[3],.combine="rbind") %dopar% (
      getcsize_base(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib)
    )
  }
  
}
