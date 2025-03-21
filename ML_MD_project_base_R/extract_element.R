extract_element <- function(elety) {
  # Match specific prefixes and patterns dynamically for elements
  if (grepl("^[A-Z]*N", elety)) return("N")  # Matches variations ending with N
  if (grepl("^[A-Z]*C", elety)) return("C")  # Matches variations ending with C
  if (grepl("^[A-Z]*O", elety)) return("O")  # Matches variations ending with O
  if (grepl("^[A-Z]*S", elety)) return("S")  # Matches variations ending with S
  return("C")  # Default if no match
}
