#!/bin/bash

# Exit on error
set -e

# Usage
show_help() {
    echo "Usage: $0 -c <config_file>"
    echo "Options: "
    echo "  -c, --config: Path to config file (YAML)"
    echo "  -h, --help: Prints help information"
    exit 1
}

# Function to extract value for a given key in YAML format
get_yaml_value() {
    local key=$1
    local default=$2
    # Isolate value per key and remove any leading and trailing/or whitespace characters (spaces, tabs etc)
    local value=$(grep "^${key}:" "$yaml_file" | cut -d':' -f2- | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    echo "${value:-$default}"
} 

# Export variables from YAML file using the custom get_yaml_value function
parse_yaml() {
    local yaml_file=$1
       
    export DATA_URL=$(get_yaml_value "data_url")
    export GENOME_URL=$(get_yaml_value "genome_url")
    export THREADS=$(get_yaml_value "threads" "1")  # Defaults to 1 if threads are not specified
}

# Function to download and extract data using wget and tar 
download_data() {
    local data_url=$1
    
    echo "Downloading data..."

    if ! command -v wget &> /dev/null; then
        echo "wget is required. Please install it first."
        exit 127
    fi

    wget -O data.tar.gz "$data_url"

    if ! command -v gunzip &> /dev/null; then
        echo "gunzip is required. Please install it first."
        exit 127
    fi

    gunzip data.tar.gz
    mkdir -p data
    tar -xf data.tar --strip-components=1 -C ./data
    rm data.tar
}

# Function to download and prepare genome
download_genome() {
    local genome_url=$1
    
    echo "Downloading genome..."

    if ! command -v wget &> /dev/null; then
        echo "wget is required. Please install it first."
        exit 127
    fi

    wget -O genome.fasta.gz "$genome_url"
    gunzip genome.fasta.gz
}

# Function to rename fastq files and keep only relevant information
rename_data() {
    local directory=$1
    echo "Renaming files..."

    if [[ -n "$(ls -A $directory/*.fastq)" ]]; then        
        for file in $(ls -1 $directory/*.fastq); do
            if [ -f "$file" ]; then
                new_name=$(echo "$file" | sed 's/.trim.sub//')
                mv "$file" "$new_name"
            fi
        done
    else 
        echo "No fastq files found in $directory".
        exit 1
    fi
}

# Function to run quality assessment using FastQC
run_fastqc() {
    local sample=$1
        
    echo "Running FastQC for sample $sample..."
    
    if ! command -v fastqc &> /dev/null; then
        echo "fastqc is required. Please install it first."
        exit 127
    fi

    fastqc "data/${sample}_1.fastq" 
    fastqc "data/${sample}_2.fastq"
}

# Function to build genome index using bwa-mem2
build_genome_index() {
    echo "Building genome index..."

    if ! command -v bwa-mem2 &> /dev/null; then
        echo "bwa-mem2 is required. Please install it first."
        exit 127
    fi

    bwa-mem2 index -p genome.fasta genome.fasta
}

# Function to perform mapping and sorting using bwa-mem2 and samtools
mapping_and_sorting() {
    local sample=$1
    local threads=${2:-1} # threads as an optional argument

    echo "Mapping and sorting sample $sample..."

    # Create temporary directory
    local temp_dir=$(mktemp -d)
    
    # Mapping using bwa-mem2
    if ! command -v bwa-mem2 &> /dev/null; then
        echo "bwa-mem2 is required. Please install it first."
        exit 127
    fi

    if ! command -v samtools &> /dev/null; then
        echo "samtools is required. Please install it first."
        exit 127
    fi

    bwa-mem2 mem -t ${threads} genome.fasta "data/${sample}_1.fastq" "data/${sample}_2.fastq" \
    | samtools view -b \
    | samtools sort -@ ${threads} -m 2G -o "${sample}.sorted.bam"
}

# Function to generate VCF using bcftools
generate_vcf() {
    local bam_files=$1
    local threads=${2:-1} # threads as an optional argument
    
    echo "Generating VCF file..."

    if ! command -v bcftools &> /dev/null; then
        echo "bcftools is required. Please install it first."
        exit 127
    fi

    # Variant calling using bcftools
    bcftools mpileup -Ou -f genome.fasta $bam_files \
    | bcftools call -Ou -mv -o all_samples.vcf
}

# Main workflow
main() {
    local config_file=""

    # Process arguments using getopt
    ARGS=$(getopt -o c:h --long config:,help -n "$0" -- "$@")

    # Stop if arguments were not processed correctly by getopt
    if [ $? != 0 ]; then
        echo "Failed to process arguments."
        exit 1
    fi

    # Evaluate arguments processed by getopt
    eval set -- "$ARGS"

    # Parse command line arguments using getopt
    while true; do
        case "$1" in 
            -c|--config)
                config_file="$2"
                shift 2
                ;;
            -h|--help)
                show_help
                exit 0
                ;;           
            --)
                shift
                break
                ;;
            *)
                echo "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done

    # If config file is empty print usage information
    if [ -z "$config_file" ]; then
        show_help
        exit 1
    fi
    
    # Parse config file - this will assign the appropriate values to the variables used downstream
    parse_yaml "$config_file"

    # Download and index genome
    if [ ! -f "genome.fasta" ]; then
        download_genome "$GENOME_URL"
        build_genome_index
    fi
    
    # Download and rename data to keep only sample information
    if [ ! -d "data" ]; then
        download_data "$DATA_URL"
    fi

    # Check if directory exists and has fastq files
    if [[ -d "data" && -n "$(ls -A data/*.trim.sub.fastq 2>/dev/null)" ]]; then
        rename_data "data"

        # Get sample names
        samples=($(ls data/*_1.fastq | sed 's/data\///g' | sed 's/_1.fastq//g'))
    fi

    # Run FastQC, mapping and sorting for all samples individually
    for sample in "${samples[@]}"; do

        # Run FastQC per sample only if a report does not already exist
        if ! [[ -f "{$sample}_fastqc.html" ]]; then
            run_fastqc "$sample" "$THREADS"
        fi

        # Run mapping and sorting unless the BAM file already exists
        if ! [[ -f "${sample}.sorted.bam" ]]; then
            mapping_and_sorting "$sample" "$THREADS"
        else
            echo "Error: ${sample}.sorted.bam is empty or was not created"
            exit 1
        fi
    done

    # Generate VCF if any BAM files were successfully created
    for bam_files in $(ls -1 *.sorted.bam); do
        echo "$bam_files"

        if [[ -s $bam_files ]]; then
            generate_vcf "$bam_files" "$THREADS"
        else
            echo "No sorted BAM files were generated"
            exit 1
        fi
    done

    # Could also be done using find
    # if [ $(find . -name "*.sorted.bam" | wc -l) -gt 0 ]; then
    # bam_files = $(ls -1 *.sorted.bam)
  
    #  for ( bam_file in "${bam_files[@]})"; then
    #       generate_vcf "${bam_files[@]}" "$THREADS"

    if [[ -s "all_samples.vcf" ]]; then
        echo "Workflow completed successfully!"   
    fi
}

# Run the workflow
main "$@"