import argparse
import pyvcf
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder
import re

parser = argparse.ArgumentParser(description='Perform PCA from input VCF file.')
parser.add_argument('-i', '--input_vcf', type=str, help='Input VCF file')
args = parser.parse_args()

def vcf2df(vcf_fname):
    """Convert a subsetted vcf file to pandas DataFrame"""
    vcf_reader = vcf.Reader(filename=vcf_fname)

    df = pd.DataFrame(index=vcf_reader.samples)
    for variant in vcf_reader:
        df[variant.ID] = [call.gt_type if call.gt_type is not None else 3 for call in variant.samples]

    return df

df = vcf2df(args.input_vcf)

ncols = len(df.columns)
ohe = OneHotEncoder(categories=[range(4)] * ncols, sparse=False)

X = ohe.fit_transform(df.values)

samples = df.index.tolist()

def run_PCA(X, n_components=3):
    pca = PCA(n_components=n_components).fit(X)
    X_red = pca.transform(X)
    df_red = pd.DataFrame(X_red, columns=[f'PCA component {i+1}' for i in range(n_components)], index=[re.sub(".bam", "", sample) for sample in samples])
    return df_red

df_dim = run_PCA(X)

# Add 'samples' column to df_dim
df_dim['samples'] = df_dim.index

# Generate random colors for each sample
random_colors = sns.color_palette("husl", n_colors=len(samples))

# Create a dictionary mapping samples to random colors
sample_colors = dict(zip(samples, random_colors))

# Ensure that sample_colors only contains keys present in df_dim
valid_sample_colors = {sample: sample_colors.get(sample, random_colors[i]) for i, sample in enumerate(df_dim.index)}

plt.figure(figsize=(8,8))
plt.title('PCA plot')
sns.scatterplot(x='PCA component 1', y='PCA component 2', data=df_dim, hue='samples', palette=valid_sample_colors)
plt.savefig('pca_plot.png')
