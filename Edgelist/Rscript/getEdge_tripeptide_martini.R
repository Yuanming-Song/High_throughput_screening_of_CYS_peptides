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
statelist<-c("dis"
             , "mon"
)
framei<-200
#clustersize bin width
binwidth<-1
outdir<-"/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tripeptide"
csizedir<-"/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Csizedist/"
for (packagename in requiredpackages) {
  # Load bio3d library, if it doesn't exist, install it
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}

registerDoParallel(cores = detectCores())


source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getedge_base.R")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getEdge.R")
# Improved function to extract element information from elety
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")

source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_from_edge.R")



# Set main directory
maindir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide/MARTINI22"

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
    for (state in statelist) {
      # Set directory paths
      simdir <- file.path(maindir, paste0("Tripeptide_", state, "_C", pos), seq)
      
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
      
      
      edgelist<-getEdge(dimerized)
      
      # Append sizehis to the appropriate data frame
      if (exists("edgelist")) {
        if (state == "dis") {
          save(edgelist,file=file.path(outdir,"dimer",paste0(seq,".rda")))
        } else {
          save(edgelist,file=file.path(outdir,"monomer",paste0(seq,".rda")))
        }
      }
      sizehis<-getdishis_from_edge(edgelist)
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
      rm(edgelist)
      rm(sizehis)
    }
    save(dimer_sizehis,file=file.path(csizedir,"dimer_cdist_tripeptide.rda"))
    save(monomer_sizehis,file=file.path(csizedir,"monomer_cdist_tripeptide.rda"))
  }
}


