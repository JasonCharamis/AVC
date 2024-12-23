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
    local yaml_file=$3

    # Isolate value per key and remove any leading and trailing whitespace
    local value=$(grep "^${key}:" "${yaml_file}" | cut -d':' -f2- | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    
    if [[ -z "$value" && -z "$default" ]]; then
        echo "Error: Required value for ${key} not found in config file" >&2
        exit 1
    fi
    
    echo "${value:-$default}"
} 

# Export variables from YAML file
parse_yaml() {
    local yaml_file=$1
    
    # Validate yaml file exists
    if [[ ! -f "${yaml_file}" ]]; then
        echo "Error: Config file ${yaml_file} not found"
        exit 1
    fi
       
    export DATA_URL=$(get_yaml_value "data_url" "" "${yaml_file}")
    export GENOME_URL=$(get_yaml_value "genome_url" "" "${yaml_file}")
    export THREADS=$(get_yaml_value "threads" "1" "${yaml_file}")

    # Validate required URLs
    if [[ -z "${DATA_URL}" || -z "${GENOME_URL}" ]]; then
        echo "Error: data_url and genome_url are required in config file"
        exit 1
    fi
}

# Function to download and extract data
download_data() {
    local data_url=$1
    local temp_dir=$(mktemp -d)
    trap 'rm -rf "$temp_dir"' EXIT
    
    echo "Downloading data..."

    if ! command -v wget &> /dev/null; then
        echo "wget is required. Please install it first."
        exit 127
    fi

    if ! wget -O "${temp_dir}/data.tar.gz" "${data_url}"; then
        echo "Error: Failed to download data from ${data_url}"
        exit 1
    fi

    if ! command -v gunzip &> /dev/null; then
        echo "gunzip is required. Please install it first."
        exit 127
    fi

    gunzip "${temp_dir}/data.tar.gz"
    mkdir -p data
    if ! tar -xf "${temp_dir}/data.tar" --strip-components=1 -C ./data; then
        echo "Error: Failed to extract data archive"
        exit 1
    fi
    rm -f "${temp_dir}/data.tar"
}

# Download and prepare genome for mapping
download_genome() {
    local genome_url=$1
    
    echo "Downloading genome..."

    if ! command -v wget &> /dev/null; then
        echo "wget is required. Please install it first."
        exit 127
    fi

    if ! wget -O genome.fasta.gz "${genome_url}"; then
        echo "Error: Failed to download genome from ${genome_url}"
        exit 1
    fi
    
    if ! gunzip genome.fasta.gz; then
        echo "Error: Failed to extract genome file"
        exit 1
    fi
}

# Rename fastq files
rename_data() {
    local directory=$1
    echo "Renaming files in ${directory}..."

    for file in "${directory}"/*.fastq; do
        if [[ -f "${file}" ]]; then
            new_name=$(echo "${file}" | sed 's/.trim.sub//')
            mv "${file}" "${new_name}"
        else
            echo "No fastq files found in ${directory}"
            exit 1
        fi
    done
}

# Run FastQC
run_fastqc() {
    local sample=$1
        
    echo "Running FastQC for sample ${sample}..."
    
    if ! command -v fastqc &> /dev/null; then
        echo "fastqc is required. Please install it first."
        exit 127
    fi

    if [[ ! -f "data/${sample}_1.fastq" || ! -f "data/${sample}_2.fastq" ]]; then
        echo "Error: Input fastq files not found for sample ${sample}"
        exit 1
    fi

    fastqc "data/${sample}_1.fastq" 
    fastqc "data/${sample}_2.fastq"
}

# Build genome index using bwa-mem2
build_genome_index() {
    echo "Building genome index..."

    if ! command -v bwa-mem2 &> /dev/null; then
        echo "bwa-mem2 is required. Please install it first."
        exit 127
    fi

    if [[ ! -f genome.fasta ]]; then
        echo "Error: genome.fasta file not found"
        exit 1
    fi

    if ! bwa-mem2 index -p genome.fasta genome.fasta; then
        echo "Error: Failed to build genome index"
        exit 1
    fi
}

# Map reads to reference genome
map_reads() {
    local sample=$1
    local threads=${2:-1}
    local output="${sample}.raw.bam"

    echo "Mapping reads for sample ${sample}..."

    if ! bwa-mem2 mem -t "${threads}" genome.fasta \
        "data/${sample}_1.fastq" "data/${sample}_2.fastq" \
        | samtools view -b > "${output}"; then
        echo "Error: Mapping failed for sample ${sample}"
        exit 1
    fi

    if [[ ! -s "${output}" ]]; then
        echo "Error: Output BAM file is empty for sample ${sample}"
        exit 1
    fi
}

# Process BAM file through collate, fixmate, sort, and markdup
process_bam() {
    local sample=$1
    local threads=${2:-1}
    local input="${sample}.raw.bam"
    
    echo "Processing BAM file for sample ${sample}..."

    # Check required tools
    if ! command -v samtools &> /dev/null; then
        echo "samtools is required. Please install it first."
        exit 127
    fi

    # Collate reads by name
    echo "Collating reads by name..."
    if ! samtools collate -@ "${threads}" "${input}" "${sample}.collate"; then
        echo "Error: Collate failed for sample ${sample}"
        exit 1
    fi

    # Add mate information
    echo "Adding mate information..."
    if ! samtools fixmate -m "${sample}.collate.bam" "${sample}.fixmate.bam"; then
        echo "Error: Fixmate failed for sample ${sample}"
        exit 1
    fi

    # Sort by coordinate
    echo "Sorting by coordinate..."
    if ! samtools sort -@ "${threads}" -m 2G "${sample}.fixmate.bam" > "${sample}.sorted.bam"; then
        echo "Error: Sort failed for sample ${sample}"
        exit 1
    fi

    # Mark and remove duplicates
    echo "Marking and removing duplicates..."
    if ! samtools markdup -r -s -@ "${threads}" "${sample}.sorted.bam" "${sample}.markdup.bam" \
        2> "${sample}.duplicate_metrics.txt"; then
        echo "Error: Markdup failed for sample ${sample}"
        exit 1
    fi

    # Clean up intermediate files
    rm -f "${sample}.collate.bam" "${sample}.fixmate.bam" "${sample}.sorted.bam"
}

# Generate VCF using bcftools
generate_vcf() {
    local -a bam_files=("${@:1:$#-2}")
    local genome="${@: -2:1}"
    local threads="${!#}"

    echo "Generating VCF..."
    echo "BAM files: ${bam_files[@]}"
    echo "Genome: ${genome}"
    echo "Threads: ${threads}"
    
    # Validate inputs
    if [[ ! -f "${genome}" ]]; then
        echo "Error: Genome file ${genome} not found"
        exit 1
    fi

    for bam in "${bam_files[@]}"; do
        if [[ ! -f "${bam}" ]]; then
            echo "Error: BAM file ${bam} not found"
            exit 1
        fi
    done
    
    if ! bcftools mpileup \
        -Ou \
        -p \
        -f "$genome" \
        --threads "$threads" \
        "${bam_files[@]}" \
        | bcftools call \
        -mv \
        -Oz \
        --threads "$threads" \
        -o all_samples.vcf; then
        echo "Error: VCF generation failed"
        exit 1
    fi

    if [[ ! -s all_samples.vcf ]]; then
        echo "Error: Output VCF file is empty or was not created"
        exit 1
    fi
}

# Run PCA
run_pca() {
    local sample_vcf=$1

    if ! command -v python &> /dev/null; then
        echo "Python is required. Please install it first."
        exit 127
    fi

    if [[ ! -f "${sample_vcf}" ]]; then
        echo "Error: VCF file ${sample_vcf} not found"
        exit 1
    fi

    echo "Running PCA on VCF file..."
    if ! python scripts/pca.py -vcf "${sample_vcf}"; then
        echo "Error: PCA analysis failed"
        exit 1
    fi
}

# Main workflow
main() {
    local config_file=""

    # Import arguments using getopt
    ARGS=$(getopt -o c:h --long config:,help -n "$0" -- "$@")

    if [[ $? != 0 ]]; then
        echo "Failed to process arguments."
        exit 1
    fi

    eval set -- "$ARGS"

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

    if [[ -z "${config_file}" ]] || [[ $# -gt 0 ]]; then
        show_help
        exit 1
    fi
    
    # Parse config file
    parse_yaml "${config_file}"

    # Download and index genome
    if [[ ! -f "genome.fasta" ]]; then
        download_genome "${GENOME_URL}"
        build_genome_index
    fi
    
    # Download and rename data
    if [[ ! -d "data" ]]; then
        download_data "${DATA_URL}"
    fi

    # Check for input fastq files
    if ls data/*.fastq >/dev/null 2>&1; then
        if [[ $(ls data/*.fastq | grep 'trim.sub' ) ]]; then
            rename_data "data"
        fi

        samples=($(ls data/*_1.fastq | sed 's/data\///g' | sed 's/_1.fastq//g'))
    else
        echo "No input fastq files found"
        exit 1
    fi

    # Validate samples array
    if [[ ${#samples[@]} -eq 0 ]]; then
        echo "Error: No valid samples found"
        exit 1
    fi

    # Process each sample
    for sample in "${samples[@]}"; do

        # Run FastQC if needed
        if [[ ! -f "data/${sample}_1_fastqc.html" ]]; then
            run_fastqc "$sample"
        fi

        # Map reads if needed
        if [[ ! -s "${sample}.raw.bam" ]]; then
            map_reads "$sample" "$THREADS"
        fi

        # Process BAM if needed
        if [[ ! -s "${sample}.markdup.bam" ]]; then
            process_bam "$sample" "$THREADS"
        fi
    done

    # Collect final BAM files
    bam_files=()
    for bamfile in *.markdup.bam; do
        if [[ -f "$bamfile" && -s "$bamfile" ]]; then
            bam_files+=("$bamfile")
        fi
    done
    
    # Generate VCF
    if [[ ! -f "all_samples.vcf" && ${#bam_files[@]} -gt 0 ]]; then
        generate_vcf "${bam_files[@]}" "genome.fasta" "$THREADS"
    elif [[ ${#bam_files[@]} -eq 0 ]]; then
        echo "Error: No valid BAM files found for VCF generation"
        exit 1
    fi

    # Run PCA
    if [[ ${#bam_files[@]} -ge 2 ]]; then
        echo "These BAM files will be used for PCA: ${bam_files[@]}"
        if [[ -f "all_samples.vcf" ]]; then
            run_pca "all_samples.vcf"
            echo "Workflow completed successfully!" 
        fi
    else
        echo "Error: Not enough input BAM files (minimum 2) for performing PCA"
        exit 1
    fi
}

# Run the workflow
main "$@"