#!/bin/bash
# Define source and destination directories
src="/dfs9/tw/yuanmis1/Rscript/ML_MD_project_base"
dest="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/ML_MD_project_base_R"

# Check if source exists
if [ ! -d "$src" ]; then
  echo "Source directory $src does not exist."
  exit 1
fi

# Check if destination already exists; exit to avoid overwriting if it does
if [ -d "$dest" ]; then
  echo "Destination directory $dest already exists. Exiting."
  exit 1
fi

# Move the directory from the source to the destination
mv "$src" "$dest"

# Create a symbolic link at the original location pointing to the new destination
ln -s "$dest" "$src"

echo "Moved $src to $dest and created symbolic link at $src."