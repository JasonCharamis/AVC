import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import gzip

def parse_vcf(vcf_path):
    """Parse VCF file and return genotype matrix."""
    genotypes = []
    sample_ids = []
    
    # Handle both gzipped and regular VCF files
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'
    
    with open_func(vcf_path, mode) as f:
        for line in f:
            # Skip comment lines
            if line.startswith('##'):
                continue
            
            # Get sample IDs from header line
            if line.startswith('#CHROM'):
                sample_ids = line.replace('.sorted.bam', '').strip().split('\t')[9:] # Clean sample IDs
                continue
            
            # Process variant lines
            fields = line.strip().split('\t')
            if len(fields) < 10:  # Skip lines without genotypes
                continue
            
            # Extract genotype values for each sample
            gts = []
            for sample in fields[9:]:
                gt = sample.split(':')[0]  # Get GT field
                if gt in ['0/0', '0|0']:
                    gts.append(0)
                elif gt in ['0/1', '1/0', '0|1', '1|0']:
                    gts.append(1)
                elif gt in ['1/1', '1|1']:
                    gts.append(2)
                else:
                    gts.append(np.nan)  # Missing data
            
            # Only keep variants with less than 10% missing data
            missing_rate = np.sum(np.isnan(gts)) / len(gts)
            if missing_rate < 0.1:
                genotypes.append(gts)
    
    return np.array(genotypes).T, sample_ids  # Transpose to get samples x variants

def perform_pca(genotype_matrix, n_components=2):
    """Perform PCA on genotype matrix."""
    # Handle missing values by replacing with mean
    col_mean = np.nanmean(genotype_matrix, axis=0)
    inds = np.where(np.isnan(genotype_matrix))
    genotype_matrix[inds] = np.take(col_mean, inds[1])
    
    # Center the data
    genotype_matrix = genotype_matrix - np.mean(genotype_matrix, axis=0)
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    pc_coordinates = pca.fit_transform(genotype_matrix)
    explained_variance = pca.explained_variance_ratio_
    
    return pc_coordinates, explained_variance

def main(vcf_path):
    """Main function to process VCF and generate PCA results."""
    print("Reading VCF file...")
    genotype_matrix, sample_ids = parse_vcf(vcf_path)
    
    print("Performing PCA...")
    pc_coords, exp_var = perform_pca(genotype_matrix)
    
    # Save results
    print("Saving results...")
    with open("pca_coordinates.txt", 'w') as f:
        f.write("Sample\tPC1\tPC2\n")
        for sample, coord in zip(sample_ids, pc_coords):
            f.write(f"{sample}\t{coord[0]}\t{coord[1]}\n")
    
    print("\nVariance explained:")
    print(f"PC1: {exp_var[0]:.2%}")
    print(f"PC2: {exp_var[1]:.2%}")
    
    # Create PCA plot
    plt.figure(figsize=(10, 8))
    plt.scatter(pc_coords[:, 0], pc_coords[:, 1], c=range(len(sample_ids)), cmap='Accent')
    
    # Add sample labels if fewer than 50 samples
    if len(sample_ids) < 50:
        for i, sample in enumerate(sample_ids):
            plt.annotate(sample, (pc_coords[i, 0], pc_coords[i, 1]), 
                        xytext=(5, 5), textcoords='offset points')

    plt.xlabel(f'PC1 ({exp_var[0]:.2%} variance explained)')
    plt.ylabel(f'PC2 ({exp_var[1]:.2%} variance explained)')
    plt.title('PCA based on VCF Genotypes',fontsize=12)
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Save plot
    plt.savefig("pca_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate PCA from VCF file")
    parser.add_argument("-vcf", "--vcf_file", type=str, help="Path to input VCF file")
    args = parser.parse_args()
    
    main(args.vcf_file)