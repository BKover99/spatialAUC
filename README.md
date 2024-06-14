# Spatial AUC
Spatial AUC is a Python package that calculates spatial autocorrelation of gene sets using Moran's I statistic. It provides a convenient way to analyze the spatial patterns of gene expression in spatial transcriptomics data.

Features

Retrieve gene sets from MSigDB database
Calculate area under the curve (AUC) for gene sets
Compute spatial autocorrelation using Moran's I statistic
Utilize Numba for performance optimization
Seamlessly integrate with Scanpy and Squidpy for spatial analysis

Installation
To install Spatial AUC, you can use pip:
shellCopy codepip install spatial-auc

Make sure you have the following dependencies installed:

numpy
pandas
scanpy
squidpy
numba
gseapy

Usage
Here's a basic example of how to use Spatial AUC:
pythonCopy code

import scanpy as sc
from spatial_auc import spatial_auc

# Load your Scanpy AnnData object
adata = sc.read_h5ad("your_data.h5ad")

# Calculate spatial autocorrelation
morans_table, adata_updated = spatial_auc(adata, gene_sets=["m5.all", "m2.all"], n_perms=1000)

# Access the results
print(morans_table)

For more detailed examples and usage instructions, please refer to the Google Colab examples.

Your dataset is too large for computations like this? Consider using [Pseudovisium](https://github.com/BKover99/Pseudovisium) for hexagonal-binning and make your analysis an order of magnitude faster and more memory-efficient.


### Author
Bence Kover
https://twitter.com/kover_bence 
https://www.linkedin.com/in/ben-kover/
