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
  ,"tidyr"
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

source("/dfs9/tw/yuanmis1/Rscript/fromAlfredo/shape.R")
logfile="csize_anal.log"

source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/gyrT_largest.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getGyrT.R")


binwidth<-1

# Read trajectory files
# Improved function to extract element information from elety




# Set main directory
maindir <- "/share/crsp/lab/dtobias/yuanmis1/mrsec/ML-MD-Peptide/MARTINI22/"

# Initialize result storage for dimer and monomer sizehis
dimer_sizehis<-seq(binwidth,300+binwidth,binwidth)
monomer_sizehis <- seq(binwidth,300+binwidth,binwidth)
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y","E","D") # Exclude "C"



gyrT<-list()
clussize<-list()
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
    gyrT[[seq]]<-list()
    clussize[[seq]]<-list()
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
      
      
      out<-getGyrT(dimerized)
      clussize[[seq]][[state]]<-c()
      GyrTtemp<-c()
      for (frame in 1:length(out)) {
        clussize[[seq]][[state]]<-rbind(clussize[[seq]][[state]],out[[frame]][[1]])
        GyrTtemp<-rbind( GyrTtemp,out[[frame]][[2]])
      }
      if (dimerized) {
        clussize[[seq]][[state]][,2]<-2*clussize[[seq]][[state]][,2]
      }
      # Columns corresponding to T1, T2, and T3.
      df <- data.frame(
        step = 1:nrow(GyrTtemp),
        T1 = GyrTtemp[,1],
        T2 = GyrTtemp[,2],
        T3 = GyrTtemp[,3]
      )
      
      # Reshape the data from wide to long format
      df_long <- pivot_longer(df, cols = c("T1", "T2", "T3"), names_to = "Series", values_to = "Value")
      
      
      # Append sizehis to the appropriate data frame
      if (exists("df_long")) {
        gyrT[[seq]][[state]]<-df_long

        rm(df_long)
      }
    }
  }
  save(clussize,file="clussize_tripeptide.rda")
  save(gyrT,file="gyrT_tripeptide.rda")
}




