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
)
framei<-200
for (packagename in requiredpackages) {
  # Load bio3d library, if it doesn't exist, install it
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}

registerDoParallel(cores = detectCores())


logfile="csize_anal.log"


getdishis_base<-function(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib) {
  #generate edge
  edges<-find_contacts_for_residue(AtomIndexLib, simtraj$coord[,,frame], cutoff_matrix, AtomIndexList, AtomTypeLib,c(max(simtraj$coord[,1,frame]), max(simtraj$coord[,2,frame]), max(simtraj$coord[,2,frame])))
  #format it properly
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  #do sna analysis
  compdist<-component.dist(edges)
  #return histogram
  table(cut(compdist$csize,breaks=seq(0,300+binwidth,binwidth)))
}



binwidth<-1
atom_types <- c("N", "C", "O", "S")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
cutoff_matrix <- matrix(.55, nrow = 4, ncol = 4, dimnames = list(atom_types, atom_types))
# Read trajectory files
getdishis<-function(dimerized) {
  
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
    AtomIndexListPerRes$Peptide$BB <- c(AtomIndexListPerRes$Peptide$BB,which(res_atoms$elety == "BB"))
    
    
    # Sidechains: All other atoms in the residue
    sidechain_atoms <- which(!res_atoms$elety %in% c("BB") )
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
  #get cluster size histogram for all frame, and add them up, return the result directly
  foreach(frame=framei:dim(simtraj$coord)[3],.combine="+") %dopar% (
    getdishis_base(frame,AtomIndexLib, cutoff_matrix, AtomIndexList, AtomTypeLib)
  )
  
}
# Improved function to extract element information from elety
extract_element <- function(elety) {
  # Match specific prefixes and patterns dynamically for elements
  if (grepl("^[A-Z]*N", elety)) return("N")  # Matches variations ending with N
  if (grepl("^[A-Z]*C", elety)) return("C")  # Matches variations ending with C
  if (grepl("^[A-Z]*O", elety)) return("O")  # Matches variations ending with O
  if (grepl("^[A-Z]*S", elety)) return("S")  # Matches variations ending with S
  return("C")  # Default if no match
}

cppFunction({'
    NumericMatrix find_contacts_for_residue(  
                                            List AtomIndexLib,  
                                            NumericMatrix particle,
                                            NumericMatrix cutoff_matrix,  
                                            List AtomIndexList,  
                                            List AtomTypeLib,
                                            NumericVector box_lengths) { // Added box_lengths input (Line 7)
                                            
      int max_resno = AtomIndexLib.size(); // set max num of res
      
      // Initialize an empty edge matrix with two columns
      NumericMatrix edge_matrix(0, 2);
      
      // Map atom types to matrix indices
      std::map<std::string, int> type_to_index = {{"N", 0}, {"C", 1}, {"O", 2}, {"S", 3}};
      
      // Loop over each subsequent residue for resno_i
      for (int resno_i = 1; resno_i <= max_resno; ++resno_i) {
        
        IntegerVector atoms_i_list = AtomIndexLib[resno_i - 1];
        IntegerVector atoms_i = Rcpp::as<IntegerVector>(Rcpp::wrap(atoms_i_list)); // Flattened atoms_i list
        
        // Loop over each subsequent residue for resno_j
        for (int resno_j = resno_i + 1; resno_j <= max_resno; ++resno_j) {
          
          IntegerVector atoms_j_list = AtomIndexLib[resno_j - 1];
          IntegerVector atoms_j = Rcpp::as<IntegerVector>(Rcpp::wrap(atoms_j_list)); // Flattened atoms_j list
          
          // Loop over each atom in atoms_i
          bool contact_found = false;
          for (int a = 0; a < atoms_i.size(); ++a) {
            int atom_i = atoms_i[a] - 1; // Convert to 0-based index
            
            if (atom_i < 0 || atom_i >= particle.nrow()) continue; // Check bounds
            
            std::string type_i;
            CharacterVector atom_names = AtomIndexList.names();
            for (int k = 0; k < atom_names.size(); ++k) {
              std::string name = Rcpp::as<std::string>(atom_names[k]);
              IntegerVector atom_list = AtomIndexList[name];
              if (std::find(atom_list.begin(), atom_list.end(), atoms_i[a]) != atom_list.end()) {
                type_i = Rcpp::as<std::string>(AtomTypeLib[name]);
                break;
              }
            }
            if (type_to_index.find(type_i) == type_to_index.end()) continue; // Convert to index
            int type_i_index = type_to_index[type_i];
            
            // Loop over each atom in atoms_j
            for (int b = 0; b < atoms_j.size(); ++b) {
              int atom_j = atoms_j[b] - 1; // Convert to 0-based index
              
              if (atom_j < 0 || atom_j >= particle.nrow()) continue; // Check bounds
              
              std::string type_j;
              for (int k = 0; k < atom_names.size(); ++k) {
                std::string name = Rcpp::as<std::string>(atom_names[k]);
                IntegerVector atom_list = AtomIndexList[name];
                if (std::find(atom_list.begin(), atom_list.end(), atoms_j[b]) != atom_list.end()) {
                  type_j = Rcpp::as<std::string>(AtomTypeLib[name]);
                  break;
                }
              }
              if (type_to_index.find(type_j) == type_to_index.end()) continue;
              int type_j_index = type_to_index[type_j];
              
              double cutoff = cutoff_matrix(type_i_index, type_j_index);
              
              // Apply PBC and compute corrected distance
              double dist2 = 0.0; // Squared distance
              for (int d = 0; d < 3; ++d) {
                double delta = particle(atom_i, d) - particle(atom_j, d); // Difference in coordinate
                // Apply minimum image convention (PBC correction) --> Added here
                delta -= round(delta / box_lengths[d]) * box_lengths[d];
                dist2 += delta * delta;
              }
              double dist = sqrt(dist2); // Corrected distance
              
              // Check if within cutoff
              if (dist <= cutoff) {
                // Expand edge_matrix and add new contact
                NumericMatrix new_edge_matrix(edge_matrix.nrow() + 1, 2);
                for (int row = 0; row < edge_matrix.nrow(); ++row) {
                  new_edge_matrix(row, 0) = edge_matrix(row, 0);
                  new_edge_matrix(row, 1) = edge_matrix(row, 1);
                }
                new_edge_matrix(edge_matrix.nrow(), 0) = resno_i;
                new_edge_matrix(edge_matrix.nrow(), 1) = resno_j;
                edge_matrix = new_edge_matrix;
                
                contact_found = true;
                break; // Stop inner loop
              }
            }
            if (contact_found) break;
          }
        }
      }
      return edge_matrix;
    }'})



# Set main directory
maindir <- "/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/MARTINI22/"

# Initialize result storage for dimer and monomer sizehis
dimer_sizehis<-seq(binwidth,300+binwidth,binwidth)
monomer_sizehis <- seq(binwidth,300+binwidth,binwidth)
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y","E","D") # Exclude "C"

# Loop over dipeptide positions
positions <- c(1, 2,3) # C can be at position 1 or 2
for (pos in positions) {
    if (pos == 1) {
      sequences <- outer(residues, residues, function(x, y) paste0("C_", x, "_", y))
  } else if (pos == 2) {
    sequences <- outer(residues, residues, function(x, y) paste0(x, "_C_", y))
  } else {
    sequences <- outer(residues, residues, function(x, y) paste0(x, "_", y, "_C"))
  }
  
  # Loop through each sequence
  for (seq in sequences) {
    for (state in c("dis", "mon")) {
      # Set directory paths
      simdir <- ifelse(pos==3 & state =="mon",file.path(maindir, paste0("Tripeptide_", state, "_C", pos,"_redo"), seq), file.path(maindir, paste0("Tripeptide_", state, "_C", pos), seq))
      
      # Find simulation files
      gro_file <- Sys.glob(file.path(simdir, "*md*gro"))
      xtc_file <- Sys.glob(file.path(simdir, "*md*xtc"))
      
      # Check if files exist
      if (length(gro_file) == 0 || length(xtc_file) == 0) {
        # Log missing files
        cat(paste0(seq, " (", state, ") Fail\n"), file = logfile, append = TRUE)
        next
      }
      
      # Load trajectory
      simtraj <- tryCatch({
        read.trj(xtc_file[1], gro_file[1])
      }, error = function(e) {
        # Log failure if loading fails
        cat(paste0(seq, " (", state, ") Fail\n"), file = logfile, append = TRUE)
        return(NULL)
      })
      
      # Check for incomplete trajectories
      if (!is.null(simtraj) && dim(simtraj$coord)[3] < 400) {
        cat(paste0(seq, " (", state, ") Incomplete\n"), file = logfile, append = TRUE)
        next
      }
      
      # Set dimerization state
      dimerized <- ifelse(state == "dis", 1, 0)
      
      
      sizehis<-getdishis(dimerized)
      sizehis<-sizehis/sum(sizehis)
      
      # Append sizehis to the appropriate data frame
      if (exists("sizehis")) {
        if (state == "dis") {
          dimer_sizehis <- cbind(dimer_sizehis, sizehis)
          colnames(dimer_sizehis)[ncol(dimer_sizehis)]<-seq
        } else {
          monomer_sizehis <- cbind(monomer_sizehis, sizehis)
          colnames(monomer_sizehis)[ncol(monomer_sizehis)]<-seq
        }

      }
    }
  }
  save(dimer_sizehis,file="dimer_cdist_tripeptide.rda")
  save(monomer_sizehis,file="monomer_cdist_tripeptide.rda")
}




