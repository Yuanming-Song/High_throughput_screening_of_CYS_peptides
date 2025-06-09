#!/bin/bash

# Define the main directory
maindir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI22"

# Function to submit job with appropriate parameters
submit_job() {
    local system=$1  # mon or dis
    local cindex_dir=$2
    local peptide_dir=$3
    local aa1=$4
    local aa2=$5
    local aa3=$6
    local aa4=$7
    local script_name=$8
    
    # Set runtime
    local runtime="24:00:00"
    
    echo -e "\nJob Details:"
    echo "Full peptide directory: ${cindex_dir}/${peptide_dir}"
    echo "Submission directory: ${cindex_dir}"
    echo "Script to be used: ${script_name}"
    echo -e "Press Enter to confirm job submission and directory removal, any other key to skip: \c"
    read -r response
    
    if [[ -z "$response" ]]; then
        echo "Removing old directory..."
        rm -rf "${cindex_dir}/${peptide_dir}"
        
        cd "$cindex_dir" || exit
        
        echo "Submitting job..."
        # Create exclude nodes string
        local exclude_nodes=$(tr '\n' ',' < "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/nodes_to_avoid_gromacs" | sed 's/,$//')
        
        sbatch --job-name="${peptide_dir}" \
            --mail-user=yuanmis1@uci.edu --mail-type=FAIL \
            --account=dtobias_lab \
            --partition=standard \
            --nodes=1 \
            --ntasks-per-node=1 \
            --cpus-per-task=40 \
            --time=$runtime \
            --exclude="$exclude_nodes" \
            --out=out/${peptide_dir}.out \
            --wrap="bash //dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Tetrapeptide/MARTINI_Setup_script/${script_name} ${aa1} ${aa2} ${aa3} ${aa4}"
        
        echo "Resubmitted job for ${peptide_dir} (${system})"
    else
        echo "Skipping ${peptide_dir}"
    fi
}

# Function to process a system type
process_system() {
    local system=$1  # mon or dis
    local cpos=$2
    
    local cindex_dir="$maindir/Tetrapeptide_${system}_C${cpos}"
    
    # Set script name based on system type
    local script_name
    if [[ "$system" == "mon" ]]; then
        script_name="Cys_unstapled_tetrapeptide_MARTINI_run.sh"
    else
        script_name="Cys_stapled_tetrapeptide_MARTINI_run.sh"
    fi
    
    # Check if directory exists
    if [[ ! -d "$cindex_dir" ]]; then
        echo "Directory $cindex_dir does not exist"
        return
    fi
    
    # Create temporary files for directory listings
    local temp_dir=$(mktemp)
    local temp_top=$(mktemp)
    local temp_gro=$(mktemp)
    
    # Get lists of directories and files
    cd "$cindex_dir" || exit
    ls -d */ 2>/dev/null | sed 's#/$##' > "$temp_dir"
    find . -maxdepth 2 -name "*.top" | cut -d'/' -f2 | sort | uniq > "$temp_top"
    find . -maxdepth 2 -name "*md*.gro" | cut -d'/' -f2 | sort | uniq > "$temp_gro"
    
    # Count files
    local total_dirs=$(wc -l < "$temp_dir")
    local top_count=$(wc -l < "$temp_top")
    local gro_count=$(wc -l < "$temp_gro")
    local failed_count=0
    
    echo -e "\nSystem: ${system}, Position: C${cpos}"
    echo "Total directories: ${total_dirs}"
    echo "Directories with .top files: ${top_count}"
    echo "Directories with .gro files: ${gro_count}"
    
    # Create array of failed jobs first
    declare -a failed_jobs
    while IFS= read -r dir; do
        if grep -q "^${dir}$" "$temp_top" && ! grep -q "^${dir}$" "$temp_gro"; then
            failed_jobs+=("$dir")
            ((failed_count++))
        fi
    done < "$temp_dir"
    
    echo "Number of failed jobs: ${failed_count}"
    
    if [ $failed_count -eq 0 ]; then
        echo "No failed jobs found for ${system} C${cpos}"
        rm -f "$temp_dir" "$temp_top" "$temp_gro"
        return
    fi
    
    # Process failed jobs
    echo -e "\nWould you like to process failed jobs for ${system} C${cpos}? (y/n): "
    read -r process_response
    
    if [[ "$process_response" != "y" ]]; then
        echo "Skipping all jobs for ${system} C${cpos}"
        rm -f "$temp_dir" "$temp_top" "$temp_gro"
        return
    fi
    
    local job_counter=1
    for dir in "${failed_jobs[@]}"; do
        echo -e "\nProcessing job ${job_counter}/${failed_count} for ${system} C${cpos}"
        echo "Found failed ${system} job in ${dir}"
        
        # Extract amino acid sequence from directory name
        IFS='_' read -r aa1 aa2 aa3 aa4 <<< "$dir"
        
        submit_job "$system" "$cindex_dir" "$dir" "$aa1" "$aa2" "$aa3" "$aa4" "$script_name"
        
        ((job_counter++))
    done
    
    # Cleanup temporary files
    rm -f "$temp_dir" "$temp_top" "$temp_gro"
}

# Process monomer system 
echo "Processing monomer system..."
#process_system "mon" "4"

# Process dimer system (C1-C4)
echo -e "\nProcessing dimer system..."
for cpos in {1..4}; do
    process_system "dis" "$cpos"
done 