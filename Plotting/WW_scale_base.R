# Define the Wimley–White octanol scale (ΔG in kcal/mol)
ww_octanol_scale <- c(
  A = 0.50, R = 1.81, N = 0.85, D = 3.64, C = -0.02,
  Q = 0.77, E = 3.63, G = 1.15, H = 2.33, I = -1.12,
  L = -1.25, K = 2.80, M = -0.67, F = -1.71, P = 0.14,
  S = 0.46, T = 0.25, W = -2.09, Y = -0.71, V = -0.46
)

# Function to calculate the total hydrophobicity of a peptide sequence
calculate_ww_hydrophobicity <- function(sequence) {
  # Convert sequence to uppercase and split into individual residues
  residues <- strsplit(toupper(sequence), "")[[1]]
  
  # Check for invalid residues
  invalid_residues <- setdiff(residues, names(ww_octanol_scale))
  if (length(invalid_residues) > 0) {
    stop(paste("Invalid residues found:", paste(invalid_residues, collapse = ", ")))
  }
  
  # Get sequence length
  seq_length <- length(residues)
  
  # Calculate theoretical min and max for normalization
  theoretical_min <- min(ww_octanol_scale) * seq_length  # most hydrophilic possible
  theoretical_max <- max(ww_octanol_scale) * seq_length  # most hydrophobic possible
  
  # Retrieve hydrophobicity values for each residue
  hydrophobicity_values <- ww_octanol_scale[residues]
  
  # Calculate total hydrophobicity
  total_hydrophobicity <- sum(hydrophobicity_values)
  
  # Calculate normalized hydrophobicity
  normalized_hydrophobicity <- (total_hydrophobicity - theoretical_min) / (theoretical_max - theoretical_min)
  
  # Return results
  list(
    sequence = sequence,
    total_hydrophobicity = total_hydrophobicity,
    normalized_hydrophobicity = normalized_hydrophobicity,
    residue_scores = hydrophobicity_values,
    theoretical_min = theoretical_min,
    theoretical_max = theoretical_max
  )
}

# Example peptide sequence
#peptide <- "ACDEFGHIKLMNPQRSTVWY"

# Calculate hydrophobicity
#result <- calculate_ww_hydrophobicity(peptide)

# Display results
#print(result$total_hydrophobicity)      # Total hydrophobicity score
#print(result$normalized_hydrophobicity)  # Normalized hydrophobicity score
#print(result$residue_scores)            # Hydrophobicity per residue