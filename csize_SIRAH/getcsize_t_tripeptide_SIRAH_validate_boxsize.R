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
framei<-1
for (packagename in requiredpackages) {
  # Load bio3d library, if it doesn't exist, install it
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
  library(packagename, character.only = TRUE)
}

registerDoParallel(cores = detectCores())

source("/dfs9/tw/yuanmis1/Rscript/fromAlfredo/shape.R")

source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_base.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/gyrT_largest.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/find_contacts_for_residue.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/extract_element.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/MARTINIcutoff.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getdishis_SIRAH.R")
source("/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base/getcsize_t.R")

binwidth<-1

# Read trajectory files
# Improved function to extract element information from elety

# Set main directory
maindir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tripeptide/SIRAH/Validation/Box_size/"
outmaindir<- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/csize_SIRAH/"
logfile=paste0(outmaindir,"csize_anal.log")

# Initialize result storage for dimer and monomer sizehis
dimer_sizehis<-seq(binwidth,300+binwidth,binwidth)
monomer_sizehis <- seq(binwidth,300+binwidth,binwidth)

# Initialize cumulative data frame
sum_size_t <- NULL

# Get list of box directories under maindir (maindir is assumed to be defined)
box_dirs <- list.files(path = maindir, pattern = "^[0-9]+$", full.names = TRUE)

# Loop over each box directory
for (box_dir in box_dirs) {
  box_size <- basename(box_dir)  # Get box size
  
  # Get list of tripeptide directories under the current box directory
  seq_dirs <- list.files(path = box_dir, pattern = "^[A-Z]_[A-Z]_[A-Z]$", full.names = TRUE)
  
  # Loop over each tripeptide directory
  for (seq_dir in seq_dirs) {
    # Construct the full simulation directory
    simdir <- file.path(box_dir, basename(seq_dir))
    seq_name <- basename(seq_dir)
    # Find simulation files
    gro_file <- Sys.glob(file.path(simdir, "*md*gro"))
    xtc_file <- Sys.glob(file.path(simdir, "*md*xtc"))
    
    # Check if files exist
    if (length(gro_file) == 0 || length(xtc_file) == 0) {
      # Log missing files
      cat(paste(seq_dir,seq_name,"\nFail\n"), file = logfile, append = TRUE)
      next
    }
    
    # Load trajectory
    simtraj <- tryCatch({
      read.trj(xtc_file[1], gro_file[1])
    }, error = function(e) {
      # Log failure if loading fails
      cat(paste(seq_dir,seq_name, "Fail\n"), file = logfile, append = TRUE)
      #return(NULL)
    })
    
    # Check for incomplete trajectories
    if (!is.null(simtraj) && dim(simtraj$coord)[3] < 500) {
      cat(paste0(seq_dir, seq_name, "Incomplete\n"), file = logfile, append = TRUE)
      #next
    }
    
    # Set dimerization state
    dimerized <- 0
    
    
    size_t<-getdishis_SIRAH(dimerized,FALSE)
    #sizehis<-sizehis/sum(sizehis)
    # Convert matrix to data frame and add additional columns
    size_t_df <- as.data.frame(size_t)
    size_t_df$seq <- seq_name
    size_t_df$box <- box_size
    
    # Append results to the cumulative data frame
    if (is.null(sum_size_t)) {
      sum_size_t <- size_t_df
    } else {
      sum_size_t <- rbind(sum_size_t, size_t_df)
    }
    
    # Save individual result for this simulation
    out_file <- file.path(outmaindir,"SIRAH_validation_boxsize_csize_t.rda")
    save(sum_size_t, file = out_file)
  }
}







