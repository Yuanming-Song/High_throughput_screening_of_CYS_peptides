#!/bin/bash

# Load GROMACS module
module load gromacs/2024.2/gcc.11.2.0-cuda.11.7.1

# Create output directories
mkdir -p AP_out AP_log

# Check and copy vdwradii.dat if not present
if [ ! -f "vdwradii.dat" ]; then
    echo "Copying vdwradii.dat..."
    cp /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_SIRAH/vdwradii.dat .
    cp /dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_SIRAH/residuetypes.dat .
else
    echo "vdwradii.dat already exists."
fi

# Function to calculate SASA for a single frame
calculate_sasa() {
    local input_file=$1
    local output_file=$2
    local selection=$3
    local box_size=$4
    local state=$5
    local sequence=$6
    local log_file="AP_log/${sequence}_box${box_size}_${state}_ini_sasa.log"
    echo $selection | gmx sasa -f "$input_file" -s "$input_file" -o "$output_file" > "$log_file" 2>&1
    awk '/^[^@#]/{print $2}' "$output_file"
}

# Function to process trajectory
process_traj() {
    local xtc_file=$1
    local gro_file=$2
    local output_file=$3
    local selection=$4
    local box_size=$5
    local state=$6
    local sequence=$7
    local log_file="AP_log/${sequence}_box${box_size}_${state}_sasa.log"
    
    # Calculate SASA for each frame
    echo $selection | gmx sasa -f "$xtc_file" -s "$gro_file" -o "$output_file" > "$log_file" 2>&1
}

# Function to find all available sequences
find_sequences() {
    local base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Box_size"
    local sequences=()
    
    # Find all box size directories
    for box_dir in "$base_dir"/*; do
        if [ -d "$box_dir" ] && [[ $(basename "$box_dir") =~ ^[0-9]+$ ]]; then
            # Look in monomer directories
            if [ -d "$box_dir/mon" ]; then
                for seq_dir in "$box_dir/mon"/*; do
                    if [ -d "$seq_dir" ]; then
                        seq_name=$(basename "$seq_dir")
                        sequences+=("$seq_name")
                    fi
                done
            fi
        fi
    done
    
    # Remove duplicates and sort
    printf "%s\n" "${sequences[@]}" | sort -u
}

# Function to find all available box sizes for a sequence
find_box_sizes() {
    local base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Box_size"
    local sequence=$1
    local box_sizes=()
    
    for box_dir in "$base_dir"/*; do
        if [ -d "$box_dir" ] && [[ $(basename "$box_dir") =~ ^[0-9]+$ ]]; then
            # Check both monomer and dimer directories
            for state in "mon" "dis"; do
                if [ -d "$box_dir/$state/$sequence" ]; then
                    box_sizes+=("$(basename "$box_dir")")
                    break
                fi
            done
        fi
    done
    
    # Sort numerically
    printf "%s\n" "${box_sizes[@]}" | sort -n
}

# Function to find largest box size with monomer structure
find_largest_box_with_monomer() {
    local base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Box_size"
    local sequence=$1
    local largest_box=""
    
    # Get all box sizes and sort in reverse order
    local box_sizes=($(find_box_sizes "$sequence" | sort -nr))
    
    for box_size in "${box_sizes[@]}"; do
        local mon_dir="$base_dir/$box_size/mon/$sequence"
        if [ -d "$mon_dir" ]; then
            if [ -n "$(find "$mon_dir" -name "*sol*ion*gro")" ]; then
                largest_box=$box_size
                break
            fi
        fi
    done
    
    echo "$largest_box"
}

# Function to process a single sequence
process_sequence() {
    local sequence=$1
    local base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Box_size"
    
    # Find largest box size with monomer structure
    local largest_box=$(find_largest_box_with_monomer "$sequence")
    if [ -z "$largest_box" ]; then
        echo "No suitable monomer structure found for sequence $sequence"
        return
    fi
    
    # Calculate initial SASA from largest box size monomer
    local mon_gro=$(find "$base_dir/$largest_box/mon/$sequence" -name "*sol*ion*gro" | head -n 1)
    if [ -z "$mon_gro" ]; then
        echo "Could not find monomer structure for sequence $sequence in box size $largest_box"
        return
    fi
    
    local sasa_ini=$(calculate_sasa "$mon_gro" "AP_out/${sequence}_ini.xvg" "1" "$largest_box" "mon" "$sequence")
    echo "Initial SASA for $sequence: $sasa_ini (from box size $largest_box)"
    
    # Process all box sizes
    local box_sizes=($(find_box_sizes "$sequence"))
    for box_size in "${box_sizes[@]}"; do
        # Process both monomer and dimer states
        for state in "mon" "dis"; do
            local workdir="$base_dir/$box_size/$state/$sequence"
            if [ -d "$workdir" ]; then
                local xtc_file=$(find "$workdir" -name "*md*xtc" | head -n 1)
                local gro_file=$(find "$workdir" -name "*sol*ion*gro" | head -n 1)
                
                if [ -n "$xtc_file" ] && [ -n "$gro_file" ]; then
                    local output_xvg="AP_out/${sequence}_${box_size}_${state}.xvg"
                    process_traj "$xtc_file" "$gro_file" "$output_xvg" "1" "$box_size" "$state" "$sequence"
                    
                    local dimerization_state="monomer"
                    if [ "$state" == "dis" ]; then
                        dimerization_state="dimer"
                    fi
                    
                    awk -v seq="$sequence" -v box="$box_size" -v sasa_ini="$sasa_ini" -v state="$dimerization_state" \
                        '/^[^@#]/{print $1, sasa_ini/$2, seq, box, state}' "$output_xvg" >> "AP_results.dat"
                fi
            fi
        done
    done
}

# Function to process CSH-CSSC system
process_cshcssc() {
    local base_dir="/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Dipeptide/SIRAH/Box_size"
    
    # Find all box sizes with CSH-CSSC
    local box_sizes=()
    for box_dir in "$base_dir"/*; do
        if [ -d "$box_dir" ] && [[ $(basename "$box_dir") =~ ^[0-9]+$ ]]; then
            if [ -d "$box_dir/CSH-CSSC" ]; then
                box_sizes+=("$(basename "$box_dir")")
            fi
        fi
    done
    
    # Sort numerically
    box_sizes=($(printf "%s\n" "${box_sizes[@]}" | sort -nr))
    
    # Find largest box size with CSH monomer structure
    local largest_box=""
    for box_size in "${box_sizes[@]}"; do
        local csh_dir="$base_dir/$box_size/CSH-CSSC/CSH"
        if [ -d "$csh_dir" ]; then
            if [ -n "$(find "$csh_dir/sirah_setup" -name "*sol*ion*gro")" ]; then
                largest_box=$box_size
                break
            fi
        fi
    done
    
    if [ -z "$largest_box" ]; then
        echo "No suitable CSH monomer structure found"
        return
    fi
    
    # Calculate initial SASA from largest box size CSH
    local csh_gro=$(find "$base_dir/$largest_box/CSH-CSSC/CSH/sirah_setup" -name "*sol*ion*gro" | head -n 1)
    if [ -z "$csh_gro" ]; then
        echo "Could not find CSH monomer structure in box size $largest_box"
        return
    fi
    
    local sasa_ini=$(calculate_sasa "$csh_gro" "AP_out/CSH_ini.xvg" "1" "$largest_box" "mon" "CSH")
    echo "Initial SASA for CSH: $sasa_ini (from box size $largest_box)"
    
    # Process all box sizes
    for box_size in "${box_sizes[@]}"; do
        # Process CSH
        local csh_dir="$base_dir/$box_size/CSH-CSSC/CSH"
        if [ -d "$csh_dir" ]; then
            local csh_xtc=$(find "$csh_dir/sirah_run" -name "*md*xtc" | head -n 1)
            local csh_gro=$(find "$csh_dir/sirah_setup" -name "*sol*ion*gro" | head -n 1)
            
            if [ -n "$csh_xtc" ] && [ -n "$csh_gro" ]; then
                local output_xvg="AP_out/CSH_${box_size}.xvg"
                process_traj "$csh_xtc" "$csh_gro" "$output_xvg" "1" "$box_size" "mon" "CSH"
                awk -v box="$box_size" -v sasa_ini="$sasa_ini" \
                    '/^[^@#]/{print $1, sasa_ini/$2, "CSH", box, "monomer"}' "$output_xvg" >> "AP_results.dat"
            fi
        fi
        
        # Process CSSC
        local cssc_dir="$base_dir/$box_size/CSH-CSSC/CSSC"
        if [ -d "$cssc_dir" ]; then
            local cssc_xtc=$(find "$cssc_dir/sirah_run" -name "*md*xtc" | head -n 1)
            local cssc_gro=$(find "$cssc_dir/sirah_setup" -name "*sol*ion*gro" | head -n 1)
            
            if [ -n "$cssc_xtc" ] && [ -n "$cssc_gro" ]; then
                local output_xvg="AP_out/CSSC_${box_size}.xvg"
                process_traj "$cssc_xtc" "$cssc_gro" "$output_xvg" "1" "$box_size" "dis" "CSSC"
                awk -v box="$box_size" -v sasa_ini="$sasa_ini" \
                    '/^[^@#]/{print $1, sasa_ini/$2, "CSH", box, "dimer"}' "$output_xvg" >> "AP_results.dat"
            fi
        fi
    done
}

# Main execution
echo "frame SASA_ratio sequence box_size dimerization_state" > AP_results.dat

# Process all dipeptide sequences
echo "Finding all available sequences..."
sequences=($(find_sequences))

echo "Processing dipeptide sequences..."
for sequence in "${sequences[@]}"; do
    echo "Processing sequence $sequence..."
    process_sequence "$sequence"
done

# Process CSH-CSSC system
echo "Processing CSH-CSSC system..."
process_cshcssc

echo "Analysis complete. Results saved in AP_results.dat" 