#!/usr/bin/env bash

# Default values
sample_name="sample"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -name)
            sample_name="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            echo "Usage: $0 -name <sample_name>"
            exit 1
            ;;
    esac
done

# Check if sample name is provided
if [ -z "$sample_name" ]; then
    echo "Error: Must specify sample name using -name parameter"
    echo "Usage: $0 -name <sample_name>"
    exit 1
fi

# Build filenames using sample name
input_file1="${sample_name}_1.fq.gz"
input_file2="${sample_name}_2.fq.gz"
output_file1="${sample_name}_1_clean.fq.gz"
output_file2="${sample_name}_2_clean.fq.gz"
output_dir="$sample_name"

# Execute SOAPnuke command
SOAPnuke filter -l 20 -q 0.1 -n 0.1 -1 "$input_file1" -2 "$input_file2" -C "$output_file1" -D "$output_file2" -o "$output_dir" -T 8