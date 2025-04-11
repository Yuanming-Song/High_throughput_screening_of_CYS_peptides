#!/bin/bash
output_file="performance_data.txt"
rm -f "$output_file"

echo "performance ngpu gpu version" > "$output_file"

for version_dir in $(find ./ -mindepth 1 -maxdepth 1 -type d); do
    version=$(basename "$version_dir")
    for gpu_type_dir in "$version_dir"/*; do
        gpu_type=$(basename "$gpu_type_dir")
        
        for gpu_qty_dir in "$gpu_type_dir"/*gpus; do
            gpu_qty=$(basename "$gpu_qty_dir" | sed 's/gpus//')
            
            # Check if there are any *out files in the directory
            if ls "$gpu_qty_dir"/*out 1> /dev/null 2>&1; then
                for replica_log in "$gpu_qty_dir"/*out; do
                    performance=$(grep -m 1 'Performance:' "$replica_log" | awk '{print $2}')
                    if [ -n "$performance" ]; then
                        echo "$performance $gpu_qty $gpu_type $version" >> "$output_file"
                    fi
                done
            fi
        done
    done
done
