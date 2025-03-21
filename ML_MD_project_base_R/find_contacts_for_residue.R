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
