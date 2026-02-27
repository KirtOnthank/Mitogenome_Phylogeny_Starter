#!/bin/bash

# Define directories
input_dir="cds_sequences"
output_dir="aligned"
mkdir -p $output_dir


# Align each FASTA file in the input directory and save to the output directory
for fasta_file in ${input_dir}/*.fa; do
    gene_name=$(basename ${fasta_file} .fa)
    output_file="${output_dir}/${gene_name}_aligned.fa"
    
    echo "Aligning ${gene_name}..."
    clustalo -i ${fasta_file} -o ${output_file} --force
    
    if [[ $? -ne 0 ]]; then
        echo "Failed to align ${fasta_file}"
        continue
    fi

    echo "Aligned ${fasta_file} and saved to ${output_file}"
done

echo "Finished aligning all files."
