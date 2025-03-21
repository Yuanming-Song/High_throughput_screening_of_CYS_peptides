# Function to generate tetrapeptide sequences based on pos and inpres
generate_sequences <- function(pos, inpres) {
  sequences <- c()
  if (pos == 1) {
    # If pos is 1: position1 = C, position2 = inpres, positions 3 and 4 vary
    for (r3 in residues) {
      for (r4 in residues) {
        seq <- paste("C", inpres, r3, r4, sep = "_")
        sequences <- c(sequences, seq)
      }
    }
  } else if (pos == 2) {
    # If pos is 2: position1 = inpres, position2 = C, positions 3 and 4 vary
    for (r3 in residues) {
      for (r4 in residues) {
        seq <- paste(inpres, "C", r3, r4, sep = "_")
        sequences <- c(sequences, seq)
      }
    }
  } else if (pos == 3) {
    # If pos is 3: position1 = inpres, position2 varies, position3 = C, position4 varies
    for (r2 in residues) {
      for (r4 in residues) {
        seq <- paste(inpres, r2, "C", r4, sep = "_")
        sequences <- c(sequences, seq)
      }
    }
  } else if (pos == 4) {
    # If pos is 4: position1 = inpres, positions 2 and 3 vary, position4 = C
    for (r2 in residues) {
      for (r3 in residues) {
        seq <- paste(inpres, r2, r3, "C", sep = "_")
        sequences <- c(sequences, seq)
      }
    }
  }
  return(sequences)
}
