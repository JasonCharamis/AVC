# ViS (Variant identification using Snakemake)
Fully automated RNAseq pipeline for identifying SNPs from RNAseq data using VarScan2.  
The pipeline is also available as a Docker image, to run all steps automatically in a few steps:

**1. Build Docker image**  
git clone https://github.com/JasonCharamis/ViS.git &&\  
cd ViS/workflow/ &&\  
sudo docker build -t automated_variant_identification:latest .  

**2. Run the automated Snakemake pipeline via Docker**  
sudo docker run -v $(pwd):/mnt/workdir -w /mnt/workdir automated_variant_identification:latest snakemake --cores 20 --snakefile ViS/workflow/Snakefile --use-conda --conda-frontend mamba --verbose  
# BiS
