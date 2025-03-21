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

source("/dfs9/tw/yuanmis1/Rscript/shape.R")
logfile="csize_anal.log"

source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/gyrT_largest.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")


binwidth<-1

# Read trajectory files
# Improved function to extract element information from elety




# Set main directory
maindir <- "/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/MARTINI21/"

# Initialize result storage for dimer and monomer sizehis
dimer_sizehis<-seq(binwidth,300+binwidth,binwidth)
monomer_sizehis <- seq(binwidth,300+binwidth,binwidth)

# Loop over dipeptide positions
positions <- c(1, 2) # C can be at position 1 or 2
for (pos in positions) {
  # Generate all possible dipeptide sequences
  residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y","E","D") # Exclude "C"
  sequences <- if (pos == 1) {
    paste0("C","_", residues) # C-X
  } else {
    paste0(residues,"_", "C") # X-C
  }
  
  # Loop through each sequence
  for (seq in sequences) {
    for (state in c("dis", "mon")) {
      # Set directory paths
      simdir <- file.path(maindir, paste0("Dipeptide_", state, "_C", pos), seq)
      
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
        save(dimer_sizehis,file="dimer_cdist.rda")
        save(monomer_sizehis,file="monomer_cdist.rda")
      }
    }
  }
}




