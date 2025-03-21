atom_types <- c("N", "C", "O", "S")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
cutoff_matrix <- matrix(.55, nrow = 4, ncol = 4, dimnames = list(atom_types, atom_types))