# Load required libraries
library(dplyr)
#source("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Plotting/WW_scale_base.R")
# List of sequences and their rationales
sequence_info <- list(
  list(seq = "CAS", rationale = "Done"),
  list(seq = "CGI", rationale = "Done"),
  list(seq = "CPV", rationale = "Done (+)"),
  list(seq = "TAC", rationale = "Done"),
  list(seq = "CH", rationale = "Dipeptide top Candidate"),
  list(seq = "YC", rationale = "Dipeptide top Candidate"),
  list(seq = "LGC", rationale = "Counterintuitive, dimer aggregates less?"),
  list(seq = "VVC", rationale = "Also counter intuitive, but all hydrophobic side chains"),
  list(seq = "CGAG", rationale = "highest dAP, highest Edge$_d$ in top dEdge"),
  list(seq = "VAAC", rationale = "top dAP, highest AP$_d$"),
  list(seq = "CAAA", rationale = "highest dEdge"),
  list(seq = "CDEF", rationale = "Best Aromatic containing candidate based on dEdge"),
  list(seq = "EYNC", rationale = "\\+ charged res at N-term, aromatic in the middle (Previous MARTINI based design rule)\\cite{Wang_2024,van_Teijlingen_2021,Frederix2011}"),
  list(seq = "CGQQ", rationale = "Polar side chain interactions"),
  list(seq = "KFDC", rationale = "KKFDD designed to assemble, charge balanced, but experimentally not \\cite{Batra_2022} Also negative control"),
  list(seq = "EEAC", rationale = "Negative control")
  #,
  #list(seq = "CKFD", rationale = "KKFDD designed to assemble, but experimentally not\\cite{Batra_2022}")
)

# Extract sequences for processing
sequences <- sapply(sequence_info, function(x) x$seq)

# Explicitly ensure CSH_cg and CSH_aa are included in the sequence list for benchmarking
sequences <- unique(c("CSH_cg", "CSH_aa", sequences))

# Track missing data frames (to report only once)
reported_missing <- list()

report_missing_once <- function(name) {
  if (!(name %in% reported_missing)) {
    warning(paste("Missing data frame:", name))
    reported_missing[[length(reported_missing) + 1]] <<- name
  }
}

# Function to determine peptide type
get_peptide_type <- function(seq) {
  len <- nchar(seq)
  if (len == 2) return("dipeptide")
  if (len == 3) return("tripeptide")
  if (len == 4) return("tetrapeptide")
  stop("Unsupported peptide length")
}

# Load gyrT objects for dipeptide and tripeptide once at the top
if (!exists("gyrT_dipeptide")) { load("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/gyrT/data/gyrT_dipeptide.rda"); gyrT_dipeptide <- gyrT; rm(gyrT) }
if (!exists("gyrT_tripeptide")) { load("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/gyrT/data/gyrT_tripeptide.rda"); gyrT_tripeptide <- gyrT; rm(gyrT) }

# Function to get data for a sequence
get_sequence_data <- function(seq) {
  # Special handling for CSH_cg and CSH_aa (benchmark/control)
  if (seq %in% c("CSH_cg", "CSH_aa")) {
    if (!exists("csh_cssc_results")) {
#      load("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats.rda")
      load("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/CSH-CSSC/csh_cssc_network_stats_500ns.rda")
      csh_cssc_results <<- final_results_500ns
    }
    result <- list(
      sequence = seq,
      peptide_type = "benchmark",
      dimer_ap = NA,
      monomer_ap = NA,
      dap = NA,
      dimer_edges = NA,
      monomer_edges = NA,
      dedge = NA,
      dimer_esp1 = NA,
      monomer_esp1 = NA,
      desp1 = NA,
      logP = NA,
      logP_norm = NA,
      dimer_A1A3 = "N/A",
      dimer_csize = "N/A"
    )
    # Get edge data
    mon_data <- csh_cssc_results[csh_cssc_results$seqname == seq & csh_cssc_results$state == "monomer" & csh_cssc_results$stat_name == "edges", ]
    dim_data <- csh_cssc_results[csh_cssc_results$seqname == seq & csh_cssc_results$state == "dimer" & csh_cssc_results$stat_name == "edges", ]
    result$monomer_edges <- mon_data$value
    result$dimer_edges <- dim_data$value
    result$dedge <- result$dimer_edges - result$monomer_edges

    # Get ESP1 data
    mon_esp1 <- csh_cssc_results[csh_cssc_results$seqname == seq & csh_cssc_results$state == "monomer" & csh_cssc_results$stat_name == "esp1", ]
    dim_esp1 <- csh_cssc_results[csh_cssc_results$seqname == seq & csh_cssc_results$state == "dimer" & csh_cssc_results$stat_name == "esp1", ]
    result$monomer_esp1 <- mon_esp1$value
    result$dimer_esp1 <- dim_esp1$value
    result$desp1 <- result$dimer_esp1 - result$monomer_esp1
    return(result)
  }

  # Only call get_peptide_type for real peptide sequences
  peptide_type <- get_peptide_type(seq)
  result <- list(
    sequence = seq,
    peptide_type = peptide_type,
    dimer_ap = NA,
    monomer_ap = NA,
    dap = NA,
    dimer_edges = NA,
    monomer_edges = NA,
    dedge = NA,
    dimer_esp1 = NA,
    monomer_esp1 = NA,
    desp1 = NA,
    logP = NA,
    logP_norm = NA,
    dimer_A1A3 = "N/A",
    dimer_csize = "N/A"
  )

  # Calculate WW scales
  tryCatch({
    ww_result <- calculate_ww_hydrophobicity(seq)
    result$logP <- ww_result$total_hydrophobicity
    result$logP_norm <- ww_result$normalized_hydrophobicity
  }, error = function(e) {
    warning(paste("Could not calculate WW scale for", seq))
  })

  # Get data based on peptide type
  if (peptide_type == "dipeptide") {
    if (exists("combined_df_ddedge_dipep")) {
      max_frame <- max(combined_df_ddedge_dipep$frame)
      
      # Get data for the sequence
      seq_data <- combined_df_ddedge_dipep[combined_df_ddedge_dipep$seqname == seq & 
                                            combined_df_ddedge_dipep$frame == max_frame, ]
      
      if (nrow(seq_data) > 0) {
        # Extract edge data
        edge_data <- seq_data[seq_data$stat_name == "edges", ]
        if (nrow(edge_data) > 0) {
          result$dimer_edges <- edge_data$value[edge_data$state == "dimer"]/300
          result$monomer_edges <- edge_data$value[edge_data$state == "monomer"]/300
          result$dedge <- result$dimer_edges - result$monomer_edges
          
          # Extract AP data from ratio column of edge data
          result$dimer_ap <- edge_data$ratio[edge_data$state == "dimer"]
          result$monomer_ap <- edge_data$ratio[edge_data$state == "monomer"]
          result$dap <- result$dimer_ap - result$monomer_ap
        }
        
        # Extract esp1 data
        esp1_data <- seq_data[seq_data$stat_name == "esp1", ]
        if (nrow(esp1_data) > 0) {
          result$dimer_esp1 <- esp1_data$value[esp1_data$state == "dimer"]/300
          result$monomer_esp1 <- esp1_data$value[esp1_data$state == "monomer"]/300
          result$desp1 <- result$dimer_esp1 - result$monomer_esp1
        }
      }

      # Get A1A3 and csize
      if (exists("gyrT_dipeptide") && !is.null(gyrT_dipeptide[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]])) {
        df <- gyrT_dipeptide[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]]
        last_frame <- df[df$step == max(df$step), ]
        if (nrow(last_frame) > 0) {
          Tvals <- last_frame$Value[match(c("T1", "T2", "T3"), last_frame$Series)]
          if (all(!is.na(Tvals))) {
            result$dimer_A1A3 <- sprintf("%.3f", max(Tvals) / min(Tvals))
          }
          if ("csize" %in% colnames(last_frame)) {
            result$dimer_csize <- as.character(last_frame$csize[1])
          }
        }
      }
    } else {
      report_missing_once("combined_df_ddedge_dipep")
    }
  } else if (peptide_type == "tripeptide") {
    if (exists("combined_df_ddedge_tripep")) {
      max_frame <- max(combined_df_ddedge_tripep$frame)
      
      # Get data for the sequence
      seq_data <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$seqname == seq & 
                                              combined_df_ddedge_tripep$frame == max_frame, ]
      
      if (nrow(seq_data) > 0) {
        # Extract edge data
        edge_data <- seq_data[seq_data$stat_name == "edges", ]
        if (nrow(edge_data) > 0) {
          result$dimer_edges <- edge_data$value[edge_data$state == "dimer"]/300
          result$monomer_edges <- edge_data$value[edge_data$state == "monomer"]/300
          result$dedge <- result$dimer_edges - result$monomer_edges
          
          # Extract AP data from ratio column of edge data
          result$dimer_ap <- edge_data$ratio[edge_data$state == "dimer"]
          result$monomer_ap <- edge_data$ratio[edge_data$state == "monomer"]
          result$dap <- result$dimer_ap - result$monomer_ap
        }
        
        # Extract esp1 data
        esp1_data <- seq_data[seq_data$stat_name == "esp1", ]
        if (nrow(esp1_data) > 0) {
          result$dimer_esp1 <- esp1_data$value[esp1_data$state == "dimer"]/300
          result$monomer_esp1 <- esp1_data$value[esp1_data$state == "monomer"]/300
          result$desp1 <- result$dimer_esp1 - result$monomer_esp1
        }
      }

      # Get A1A3 and csize
      if (exists("gyrT_tripeptide") && !is.null(gyrT_tripeptide[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]])) {
        df <- gyrT_tripeptide[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]]
        last_frame <- df[df$step == max(df$step), ]
        if (nrow(last_frame) > 0) {
          Tvals <- last_frame$Value[match(c("T1", "T2", "T3"), last_frame$Series)]
          if (all(!is.na(Tvals))) {
            result$dimer_A1A3 <- sprintf("%.3f", max(Tvals) / min(Tvals))
          }
          if ("csize" %in% colnames(last_frame)) {
            result$dimer_csize <- as.character(last_frame$csize[1])
          }
        }
      }
    } else {
      report_missing_once("combined_df_ddedge_tripep")
    }
  } else if (peptide_type == "tetrapeptide") {
    # Extract main stats from combined_df_netstat_tetra
    if (exists("combined_df_netstat_tetra")) {
      max_frame <- max(combined_df_netstat_tetra$frame)
      seq_data <- combined_df_netstat_tetra[combined_df_netstat_tetra$seqname == seq & combined_df_netstat_tetra$frame == max_frame, ]
      if (nrow(seq_data) > 0) {
        edge_data <- seq_data[seq_data$stat_name == "edges", ]
        if (nrow(edge_data) > 0) {
          result$dimer_edges <- edge_data$value[edge_data$state == "dimer"]/300
          result$monomer_edges <- edge_data$value[edge_data$state == "monomer"]/300
          result$dedge <- result$dimer_edges - result$monomer_edges
          result$dimer_ap <- edge_data$ratio[edge_data$state == "dimer"]
          result$monomer_ap <- edge_data$ratio[edge_data$state == "monomer"]
          result$dap <- result$dimer_ap - result$monomer_ap
        }
        esp1_data <- seq_data[seq_data$stat_name == "esp1", ]
        if (nrow(esp1_data) > 0) {
          result$dimer_esp1 <- esp1_data$value[esp1_data$state == "dimer"]/300
          result$monomer_esp1 <- esp1_data$value[esp1_data$state == "monomer"]/300
          result$desp1 <- result$dimer_esp1 - result$monomer_esp1
        }
      }
    } else {
      report_missing_once("combined_df_netstat_tetra")
    }
    # Minimal block: extract A1A3 and csize from gyrT
    seq_chars <- strsplit(seq, "")[[1]]
    pos <- which(seq_chars == "C")[1]
    inpres <- seq_chars[which(seq_chars != "C")][1]
    tetrafn <- file.path("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/gyrT/data/Tetrapeptide", sprintf("gyrT_tetrapeptide_%d_%s_single_node.rda", pos, inpres))
    if (file.exists(tetrafn)) {
      load(tetrafn)
      if (exists("gyrT") && !is.null(gyrT[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]])) {
        df <- gyrT[[paste(strsplit(seq, "")[[1]], collapse = "_")]][["dis"]]
        last_frame <- df[df$step == max(df$step), ]
        if (nrow(last_frame) > 0) {
          Tvals <- last_frame$Value[match(c("T1", "T2", "T3"), last_frame$Series)]
          if (all(!is.na(Tvals))) result$dimer_A1A3 <- sprintf("%.3f", max(Tvals) / min(Tvals))
          if ("csize" %in% colnames(last_frame)) result$dimer_csize <- as.character(last_frame$csize[1])
        }
      }
      rm(gyrT)
    }
  }
  
  return(result)
}

# Process all sequences
all_results <- lapply(sequences, get_sequence_data)

# Load CSH AP values from shell output
ap_vals <- scan("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/CSH-CSSC/CSH_AP.txt", what = list(character(), numeric(), numeric()), quiet = TRUE)
APmon <- ap_vals[[2]][1]
APdim <- ap_vals[[3]][1]

# Assign to CSH_cg (monomer and dimer AP)
for (i in seq_along(all_results)) {
  if (all_results[[i]]$sequence == "CSH_cg") {
    all_results[[i]]$monomer_ap <- APmon
    all_results[[i]]$dimer_ap <- APdim
    all_results[[i]]$dap <- APdim - APmon
  }
}

# Ensure CSH_cg and CSH_aa are first in the table
csh_cg_idx <- which(sapply(all_results, function(x) x$sequence == "CSH_cg"))
csh_aa_idx <- which(sapply(all_results, function(x) x$sequence == "CSH_aa"))
other_idx <- setdiff(seq_along(all_results), c(csh_cg_idx, csh_aa_idx))
all_results <- c(all_results[csh_cg_idx], all_results[csh_aa_idx], all_results[other_idx])

# Create LaTeX table with new format
latex_table <- c(
  "\\clearpage",
  "\\begin{sidewaystable}",
  "    \\centering",
  "    \\resizebox{\\textwidth}{!}{%",
  "    \\begin{tabular}{|c|p{4cm}|c|c|c|c|c|c|c|c|c|c|c|c|c|}",
  "        \\hline",
  # Correct column order: Sequence, Rationale, dAP, dEdge, A1/A3, Largest cluster size, dESP1, logP, logP', AP_m, AP_d, ESP1_m, ESP1_d, Edge_m, Edge_d
  "        Sequence & Rationale & \\textbf{dAP} & \\textbf{dEdge} & A1/A3 & Largest cluster size & \\textbf{dESP1} & logP & logP$^\\prime$ & AP$_m$ & AP$_d$ & ESP1$_m$ & ESP1$_d$ & Edge$_m$ & Edge$_d$  \\\\",
  "        \\hline"
)

is_italic <- FALSE
for (result in all_results) {
  # Set rationale for CSH_cg and CSH_aa
  if (result$sequence %in% c("CSH_cg", "CSH_aa")) {
    rationale <- "Previous benchmark"
    na_str <- "N/A"
  } else {
    rationale <- sapply(sequence_info, function(x) if(x$seq == result$sequence) x$rationale else NULL)
    rationale <- rationale[!sapply(rationale, is.null)][[1]]
    na_str <- "??"
  }

  # Format each value, handling NA cases
  ap_d <- ifelse(is.na(result$dimer_ap), na_str, sprintf("\\scriptsize{%.3f}", result$dimer_ap))
  ap_m <- ifelse(is.na(result$monomer_ap), na_str, sprintf("\\scriptsize{%.3f}", result$monomer_ap))
  dap <- ifelse(is.na(result$dap), na_str, sprintf("\\textbf{%.3f}", result$dap))
  edge_d <- ifelse(is.na(result$dimer_edges), na_str, sprintf("\\scriptsize{%.2f}", result$dimer_edges))
  edge_m <- ifelse(is.na(result$monomer_edges), na_str, sprintf("\\scriptsize{%.2f}", result$monomer_edges))
  dedge <- ifelse(is.na(result$dedge), na_str, sprintf("\\textbf{%.2f}", result$dedge))
  esp1_d <- ifelse(is.na(result$dimer_esp1), na_str, sprintf("\\scriptsize{%.3f}", result$dimer_esp1))
  esp1_m <- ifelse(is.na(result$monomer_esp1), na_str, sprintf("\\scriptsize{%.3f}", result$monomer_esp1))
  desp1 <- ifelse(is.na(result$desp1), na_str, sprintf("\\textbf{%.3f}", result$desp1))
  logP <- ifelse(is.na(result$logP), na_str, sprintf("\\scriptsize{%.3f}", result$logP))
  logP_norm <- ifelse(is.na(result$logP_norm), na_str, sprintf("\\scriptsize{%.3f}", result$logP_norm))

  # Use LaTeX subscript formatting for CSH_cg and CSH_aa in the sequence column, and cite Song24
  seq_display <- result$sequence
  if (seq_display == "CSH_cg") seq_display <- "CSH$_{\\mathrm{cg}}$~\\cite{Song24}"
  if (seq_display == "CSH_aa") seq_display <- "CSH$_{\\mathrm{aa}}$~\\cite{Song24}"

  # Alternate italic for sequence and rationale
  if (is_italic) {
    seq_display <- sprintf("\\textit{%s}", seq_display)
    rationale <- sprintf("\\textit{%s}", rationale)
  }
  is_italic <- !is_italic

  # New column order: Sequence, Rationale, dAP, dEdge, dESP1, logP, logP', AP_m, AP_d, ESP1_m, ESP1_d, Edge_m, Edge_d
  a1a3 <- ifelse(is.null(result$dimer_A1A3) || is.na(result$dimer_A1A3), na_str, sprintf("%s", result$dimer_A1A3))
  csize <- ifelse(is.null(result$dimer_csize) || is.na(result$dimer_csize), na_str, sprintf("%s", result$dimer_csize))
  row <- sprintf("        %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                 seq_display,
                 rationale,
                 dap,
                 dedge,
                 a1a3,
                 csize,
                 desp1,
                 logP,
                 logP_norm,
                 ap_m,
                 ap_d,
                 esp1_m,
                 esp1_d,
                 edge_m,
                 edge_d
  )
  latex_table <- c(latex_table, row)
}
print(latex_table)
# Add closing lines and update caption for aa/cg meaning and citation
latex_table <- c(latex_table,
                 "        \\hline",
                 "    \\end{tabular}%",
                 "    }",
                 "    \\caption{CSH$_{\\mathrm{cg}}$ is SIRAH force field at 50 mM, data taken from 500 ns trajectory. CSH$_{\\mathrm{aa}}$ means all-atomistic, CSH$_{\\mathrm{cg}}$ means coarse-grained \\cite{Song24}. Subscripts $m$ and $d$ indicate monomer and dimer states respectively. Edge and ESP1 values are normalized to per monomer or monomer-equivalent. logP$^\\prime$ is the normalized Wimley-White hydrophobicity scale calculated as (logP - logP$_{min}$)/(logP$_{max}$ - logP$_{min}$), where logP$_{min}$ and logP$_{max}$ are the theoretical minimum and maximum values possible for the peptide length. Bold columns (dAP, dESP1, dEdge) show differences between dimer and monomer states. \\textbf{A1/A3 and Largest cluster size are from dimer simulation.}}",
                 "    \\label{tab:exp_validation_round2}",
                 "\\end{sidewaystable}",
                 "\\clearpage")

# Save the table
if (exists("save_plots") && save_plots) {
  writeLines(latex_table,
             paste0(pltsavedir, "exp_validation_round2.tex")
  )
} 

