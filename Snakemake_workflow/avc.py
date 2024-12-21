from prefect import task, flow
import os
import re
import subprocess
from glob import glob
from typing import List, Optional
from prefect.utilities.filesystem import create_temporary_file
import yaml

@task
def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

@task
def download_data(data_url: str) -> None:
    """Download and extract data from URL."""
    subprocess.run(["wget", "-O", "data.tar.gz", data_url], check=True)
    subprocess.run(["gunzip", "data.tar.gz"], check=True)
    os.makedirs("data", exist_ok=True)
    subprocess.run(["tar", "-xf", "data.tar", "--strip-components=1", "-C", "./data"], check=True)

@task
def rename_data(files: List[str]) -> None:
    """Rename files by removing .trim.sub from filenames."""
    for file in files:
        base_name = os.path.basename(file)
        new_name = os.path.join(os.path.dirname(file), base_name.replace('.trim.sub', ''))
        os.rename(file, new_name)

@task
def check_filenames(directory_path: str, regex_pattern: str) -> Optional[List[str]]:
    """Check for files matching the regex pattern in the directory."""
    try:
        if not os.path.exists(directory_path) or not os.path.isdir(directory_path):
            print(f"The directory {directory_path} does not exist.")
            return None
        
        pattern = re.compile(regex_pattern)
        matching_filenames = [
            filename for filename in glob(os.path.join(directory_path, "*.fastq")) 
            if re.search(pattern, filename)
        ]
        
        return matching_filenames if matching_filenames else None
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

@task
def download_genome(genome_url: str) -> str:
    """Download and decompress genome file."""
    subprocess.run(["wget", "-O", "genome.fasta.gz", genome_url], check=True)
    subprocess.run(["gunzip", "genome.fasta.gz"], check=True)
    return "genome.fasta"

@task
def run_fastqc(sample: str, threads: int = 1) -> str:
    """Run FastQC on paired-end reads."""
    r1 = f"data/{sample}_1.fastq"
    r2 = f"data/{sample}_2.fastq"
    os.makedirs("fastqc", exist_ok=True)
    
    subprocess.run(["fastqc", r1, "-t", str(threads)], check=True)
    subprocess.run(["fastqc", r2, "-t", str(threads)], check=True)
    
    return f"fastqc/{sample}_fastqc/fastqc_report.html"

@task
def build_genome_index(genome_file: str) -> str:
    """Build genome index using bwa-mem2."""
    subprocess.run(["bwa-mem2", "index", "-p", genome_file, genome_file], check=True)
    with open("index_chkp", "w") as f:
        f.write("Index built successfully")
    return "index_chkp"

@task
def mapping_and_sorting(sample: str, genome_file: str, index_chkp: str) -> str:
    """Map reads to genome and sort the resulting BAM file."""
    r1 = f"data/{sample}_1.fastq"
    r2 = f"data/{sample}_2.fastq"
    output_bam = f"{sample}.bam"
    
    with create_temporary_file() as temp_bam:
        # Mapping
        bwa_cmd = f"bwa-mem2 mem {genome_file} {r1} {r2}"
        picard_convert_cmd = f"picard SamFormatConverter I=/dev/stdin O={temp_bam}"
        
        bwa_proc = subprocess.Popen(
            bwa_cmd.split(), 
            stdout=subprocess.PIPE
        )
        subprocess.run(
            picard_convert_cmd.split(), 
            stdin=bwa_proc.stdout, 
            check=True
        )
        
        # Sorting
        subprocess.run([
            "picard", "SortSam",
            f"I={temp_bam}",
            f"O={output_bam}",
            "SORT_ORDER=coordinate",
            "CREATE_INDEX=true",
            "VALIDATION_STRINGENCY=LENIENT"
        ], check=True)
    
    return output_bam

@task
def generate_vcf(bam_files: List[str], genome_file: str) -> str:
    """Generate VCF file from BAM files using bcftools."""
    output_vcf = "all_samples.vcf"
    
    mpileup_cmd = f"bcftools mpileup -Ou -f {genome_file} {' '.join(bam_files)}"
    call_cmd = f"bcftools call -Ou -mv -o {output_vcf}"
    
    mpileup_proc = subprocess.Popen(
        mpileup_cmd.split(), 
        stdout=subprocess.PIPE
    )
    subprocess.run(
        call_cmd.split(), 
        stdin=mpileup_proc.stdout, 
        check=True
    )
    
    return output_vcf

@flow
def dna_sequencing_workflow(config_path: str):
    """Main workflow for DNA sequencing analysis."""
    # Load configuration
    config = load_config(config_path)
    
    # Download and prepare data
    data_dir = "data"
    if not os.path.exists(data_dir) or not os.path.isdir(data_dir):
        download_data(config["data_url"])
    
    matching_files = check_filenames(data_dir, "trim")
    if matching_files:
        rename_data(matching_files)
    
    # Get sample names
    samples = sorted([
        f[:-8] for f in os.listdir(data_dir) 
        if f.endswith(".fastq")
    ])
    
    # Download and index genome
    genome_file = download_genome(config["genome_url"])
    index_chkp = build_genome_index(genome_file)
    
    # Run FastQC for all samples
    fastqc_results = []
    for sample in samples:
        fastqc_results.append(run_fastqc(sample, threads=config.get("threads", 1)))
    
    # Mapping and sorting for all samples
    bam_files = []
    for sample in samples:
        bam_file = mapping_and_sorting(sample, genome_file, index_chkp)
        bam_files.append(bam_file)
    
    # Generate final VCF
    final_vcf = generate_vcf(bam_files, genome_file)
    
    return final_vcf

if __name__ == "__main__":
    # Run the workflow using a YAML config file
    dna_sequencing_workflow("config.yaml")