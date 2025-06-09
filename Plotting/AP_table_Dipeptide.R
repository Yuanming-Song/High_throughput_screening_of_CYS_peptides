# Read MARTINI and SIRAH SASA result files
martini_file <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_M21/SASA_result_with_common_mon_ini.txt"
sirah_file <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_SIRAH/SASA_result_with_common_mon_ini.txt"

martini <- read.table(martini_file, header = FALSE, stringsAsFactors = FALSE)
sirah <- read.table(sirah_file, header = FALSE, stringsAsFactors = FALSE)

# Assign column names for clarity
colnames(martini) <- c("Res1", "Res2", "MAP_m", "MAP_d")
colnames(sirah)   <- c("Res1", "Res2", "SAP_m", "SAP_d")

# Get all unique non-C residues
all_residues <- unique(c(martini$Res1, martini$Res2))
X_residues <- setdiff(all_residues, "C")

# Build the table
final_table <- data.frame(
  X = X_residues,
  MAP_m_C1 = NA, MAP_d_C1 = NA, SAP_m_C1 = NA, SAP_d_C1 = NA,
  MAP_m_C2 = NA, MAP_d_C2 = NA, SAP_m_C2 = NA, SAP_d_C2 = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(X_residues)) {
  X <- X_residues[i]
  # C1: C X
  martini_C1 <- martini[martini$Res1 == "C" & martini$Res2 == X, ]
  sirah_C1   <- sirah[sirah$Res1 == "C" & sirah$Res2 == X, ]
  # C2: X C
  martini_C2 <- martini[martini$Res1 == X & martini$Res2 == "C", ]
  sirah_C2   <- sirah[sirah$Res1 == X & sirah$Res2 == "C", ]
  
  if (nrow(martini_C1) == 1) {
    final_table$MAP_m_C1[i] <- martini_C1$MAP_m
    final_table$MAP_d_C1[i] <- martini_C1$MAP_d
  }
  if (nrow(sirah_C1) == 1) {
    final_table$SAP_m_C1[i] <- sirah_C1$SAP_m
    final_table$SAP_d_C1[i] <- sirah_C1$SAP_d
  }
  if (nrow(martini_C2) == 1) {
    final_table$MAP_m_C2[i] <- martini_C2$MAP_m
    final_table$MAP_d_C2[i] <- martini_C2$MAP_d
  }
  if (nrow(sirah_C2) == 1) {
    final_table$SAP_m_C2[i] <- sirah_C2$SAP_m
    final_table$SAP_d_C2[i] <- sirah_C2$SAP_d
  }
}

# Save the table to CSV
write.csv(final_table, file = "~/Downloads/AP_table_Dipeptide_ff.csv", row.names = FALSE)

# View the table in RStudio
View(final_table) 