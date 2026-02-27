#!/bin/bash

# Define paths relative to the script's directory
SCRIPT_DIR="$(dirname "$0")"
GENOME_FASTA="$SCRIPT_DIR/../Octopus_rubescens_mitogenome.fas"
GFF_FILE="$SCRIPT_DIR/../Octopus_rubescens_mitogenome.gff"
CDS_DIR="$SCRIPT_DIR/../cds_sequences"
MITO_CDS="$SCRIPT_DIR/../mitochondrial_cds.fasta"

# Define gene name mappings (GFF name -> Standardized name)
declare -A GENE_NAME_MAP
GENE_NAME_MAP=(
    ["cob"]="cytb"
    ["nd1"]="nad1"
    ["nd2"]="nad2"
    ["nd3"]="nad3"
    ["nd4"]="nad4"
    ["nd4l"]="nad4l"
    ["nd5"]="nad5"
    ["nd6"]="nad6"
)

# Debugging: Print out the current directory and its contents
echo "Current working directory: $(pwd)"
echo "Files in current directory:"
ls -l "$SCRIPT_DIR/.."

# Check if input files exist before proceeding
if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA file not found at $GENOME_FASTA"
    exit 1
fi

if [ ! -f "$GFF_FILE" ]; then
    echo "Error: GFF file not found at $GFF_FILE"
    exit 1
fi

rm -f "$MITO_CDS"

echo "Extracting mitochondrial exon sequences from annotation file..."

# Extract exon regions from the genome FASTA using the GFF file
awk '$3 == "exon"' "$GFF_FILE" | while read -r line; do
    seq_id=$(echo "$line" | awk '{print $1}')
    start=$(echo "$line" | awk '{print $4}')
    end=$(echo "$line" | awk '{print $5}')
    strand=$(echo "$line" | awk '{print $7}')
    gene_name=$(echo "$line" | awk -F'=' '{print $NF}' | cut -d';' -f1)

    # Standardize gene names using the mapping dictionary
    if [ "${GENE_NAME_MAP[$gene_name]+_}" ]; then
        original_name=$gene_name
        gene_name="${GENE_NAME_MAP[$gene_name]}"
        echo "Mapping gene name $original_name to $gene_name"
    fi

    # Debugging: Print each exon being processed
    echo "Processing exon: $seq_id $start-$end ($strand) -> $gene_name"

    # Extract sequence using samtools
    seq=$(samtools faidx "$GENOME_FASTA" "$seq_id:$start-$end" 2>/dev/null)

    # If sequence extraction fails, print an error
    if [ -z "$seq" ]; then
        echo "Warning: No sequence found for $seq_id:$start-$end"
        continue
    fi

    # If sequence is on the reverse strand, reverse complement it
    if [ "$strand" == "-" ]; then
        seq=$(echo "$seq" | tail -n +2 | tr -d '\n' | rev | tr 'ATGCNatgcn' 'TACGNtacgn')
    else
        seq=$(echo "$seq" | tail -n +2 | tr -d '\n')
    fi

    # Format the output in FASTA format (Keep gene name in the header)
    echo ">Octopus_rubescens_$gene_name" >> "$MITO_CDS"
    echo "$seq" >> "$MITO_CDS"
done

# If no exons were extracted, print a message
if [ ! -s "$MITO_CDS" ]; then
    echo "Warning: No mitochondrial exons extracted. Check GFF file format."
    exit 1
fi

echo "Mitochondrial exon extraction complete."

# Append extracted exons to their corresponding FASTA files
echo "Appending mitochondrial exon sequences to existing FASTA files..."

# Initialize skip_append as false at the start
skip_append=false

while read -r line; do
    if [[ $line == ">"* ]]; then
        # Reset skip_append for each new sequence
        skip_append=false
        
        # Extract gene name from the header
        gene_name=$(echo "$line" | cut -d'_' -f3)

        # Standardize gene name if needed
        if [ "${GENE_NAME_MAP[$gene_name]+_}" ]; then
            gene_name="${GENE_NAME_MAP[$gene_name]}"
        fi

        fasta_file="$CDS_DIR/${gene_name}.fa"

        # Check if the file already exists
        if [ -f "$fasta_file" ]; then
            # Print which file we are appending to
            echo "Appending Octopus_rubescens to $fasta_file"

            # Write header without gene name to the appropriate file
            echo ">Octopus_rubescens" >> "$fasta_file"
        else
            # Skip appending sequence if the file doesn't exist
            echo "Skipping $gene_name: No existing file found in $CDS_DIR"
            skip_append=true
        fi
    elif [ "$skip_append" = false ]; then
        # Write sequence only if we are not skipping
        echo "$line" >> "$fasta_file"
    fi
done < "$MITO_CDS"

echo "Processing complete. Mitochondrial exons added to their corresponding FASTA files."
