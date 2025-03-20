#!/bin/bash
output_file="performance_data.txt"
rm -f "$output_file"

echo "performance option version" > "$output_file"

for version_dir in $(find ./ -mindepth 1 -maxdepth 1 -type d); do
    version=$(basename "$version_dir")
    
    for option_dir in "$version_dir"/*; do
        option=$(basename "$option_dir")

        # Check if there are any *out files in the directory
        if ls "$option_dir"/*out 1> /dev/null 2>&1; then
            for replica_log in "$option_dir"/*out; do
                performance=$(grep -m 1 'Performance:' "$replica_log" | awk '{print $2}')
                if [ -n "$performance" ]; then
                    echo "$performance $option $version" >> "$output_file"
                fi
            done
        fi
    done
done