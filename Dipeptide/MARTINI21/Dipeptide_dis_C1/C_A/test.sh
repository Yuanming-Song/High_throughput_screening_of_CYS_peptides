#!/bin/bash

# Parse options
while getopts ":b:" opt; do
  case $opt in
    b)
      base="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if the -base option is set
if [ -z "$base" ]; then
  echo "Error: -base option is required."
  exit 1
fi

# Print the base value
echo "Base: $base"
