#!/bin/bash

# Exit on error
set -e

# Usage
usage() {
    echo "Usage: $0 -c <config_file>"
    echo "  -c: Path to config file (YAML)"
    exit 1
}

# Function to parse YAML config using yq, a command-line YAML processor
parse_yaml() {
    local yaml_file=$1

    if ! command -v yq &> /dev/null; then
        echo "yq is required to parse YAML. Please install it first."
        exit 127
    fi
    
    # Export variables from YAML
    export DATA_URL=$(yq eval '.data_url' "$yaml_file")
    export GENOME_URL=$(yq eval '.genome_url' "$yaml_file")
    export THREADS=$(yq eval '.threads // "1"' "$yaml_file")
}

# Function to download and extract data using wget and tar 
download_data() {
    local data_url=$1
    
    echo "Downloading data..."

    if ! command -v wget &> /dev/null; then
        echo "wget is required. Please install it first."
        exit 127
    fi

    # Download fastq data only if they do not already exist
    if ! ls data.tar* >/dev/null 2>&1; then
        if ! command -v gunzip &> /dev/null; then
            echo "gunzip is required. Please install it first."
            exit 127
        fi
        
        wget -O data.tar.gz "$data_url"
        gunzip data.tar.gz
        mkdir -p data
    
    else
        echo "Data already downloaded."
    fi

    # Extract data from tar file only if the data directory does not exist and it's not empty
    if ! ls -d data &> /dev/null || [ -z "$(ls -A data)" ]; then
        echo "Extracting data..."
        mkdir -p data
        tar -xf data.tar --strip-components=1 -C ./data
        rm data.tar
    fi
}

# Function to rename FASTQ files
rename_data() {
    echo "Renaming files..."

    # Rename fastq files with "_R1" and "_R2" suffixes
    for file in data/*.trim.sub*.fastq; do
        if [ -f "$file" ]; then
            new_name=$(echo "$file" | sed 's/.trim.sub//')
            mv "$file" "$new_name"
        else
            echo "File $file does not exist".
        fi
    done
}


# Function to download and prepare genome
download_genome() {
    local genome_url=$1
    
    # Download genome.fasta only if it does not already exist
    if ! ls genome.fasta* >/dev/null 2>&1; then

        if ! command -v wget &> /dev/null; then
            echo "wget is required. Please install it first."
            exit 127
        fi

        echo "Downloading genome..."
        wget -O genome.fasta.gz "$genome_url"
        gunzip genome.fasta.gz
    fi

    else
        echo "Genome fasta already downloaded."
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

    if ! ls genome.fasta.amb genome.fasta.ann genome.fasta.bwt.2bit.64 genome.fasta.pac genome.fasta.0123 &> /dev/null; then
        if ! command -v bwa-mem2 &> /dev/null; then
            echo "bwa-mem2 is required. Please install it first."
            exit 127
        fi

        bwa-mem2 index -p genome.fasta genome.fasta

    else
        echo "Genome already indexed."
    fi 
}

# Function to perform mapping and sorting using bwa-mem2 and samtools
mapping_and_sorting() {
    local sample=$1
    local threads=${2:-1} # threads as an optional argument

    echo "Mapping and sorting sample $sample..."

    # Run mapping and sortings using samtools only for samples which do not have bam files
    if ! ls "${sample}.sorted.bam" &> /dev/null; then
    
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

    else 
        echo "${sample}.sorted.bam is already present. Skipping..."
    fi
}

# Function to generate VCF using bcftools
generate_vcf() {
    local bam_files=$1
    local threads=${2:-1} # threads as an optional argument

    # Variant calling using bcftools only if VCF file does not exist
    if ! ls all_samples.vcf &> /dev/null; then

        if ! command -v bcftools &> /dev/null; then
            echo "bcftools is required. Please install it first."
            exit 127
        fi

        echo "Generating VCF file..."
    
        bcftools mpileup -Ou -f genome.fasta $bam_files \
        | bcftools call -Ou -mv -o all_samples.vcf
    else 
        echo "all_samples.vcf is already present. Skipping..."
    fi
}

# Main workflow
main() {
    local config_file=""

    # Process arguments using getopt
    ARGS=$(getopt -o c: --long yaml_file: -n avc.sh -- "$@")

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
            -c|--yaml_file)
            config_file="$2"
            shift 2
            ;;
            --)
            shift
            break
            ;;
        esac
    done
    
    # If config file is empty print usage information
    if [ -z "$config_file" ]; then
        usage
    fi
    
    # Parse config file - this will assign the appropriate values to the variables used downstream
    parse_yaml "$config_file"

    # Download and index genome
    if [ ! -f "genome.fasta" ]; then
        download_genome "$GENOME_URL"
        build_genome_index
    fi
    
    # Download and prepare data if needed
    if [ ! -d "data" ]; then
        download_data "$DATA_URL"
        rename_data 
    fi
        
    # Get sample names from renamed files
    samples=($(ls data/*_1.fastq | sed 's/data\///g' | sed 's/_1.fastq//g'))
    
    # Run FastQC for all samples
    for sample in "${samples[@]}"; do
        run_fastqc "$sample" "$THREADS"
    done
    
    # Mapping and sorting for all samples
    bam_files=""

    for sample in "${samples[@]}"; do
        mapping_and_sorting "$sample" "$THREADS"
        bam_files+=" ${sample}.sorted.bam"
    done
    
    # Generate final VCF
    generate_vcf "$bam_files" "$THREADS"
    
    echo "Workflow completed successfully!"
}

# Run the workflow
main "$@"