#!/bin/bash
output_file="performance_data.txt"
rm -f $output_file

echo "performance ncore version" > "$output_file"

for version_dir in $(find ./ -mindepth 1 -maxdepth 1 -type d); do
    version=$(basename "$version_dir")
    for core_dir in "$version_dir"/*cores; do
        core_count=$(basename "$core_dir" | sed 's/cores//')

        # Check if there are any *out files in the directory
        if ls "$core_dir"/*out 1> /dev/null 2>&1; then
            for replica_log in "$core_dir"/*out; do
                performance=$(grep -m 1 'Performance:' "$replica_log" | awk '{print $2}')
                if [ -n "$performance" ]; then
                    echo "$performance $core_count $version" >> "$output_file"
                fi
            done
        fi
    done
done