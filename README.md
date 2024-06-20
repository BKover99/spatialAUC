# spatialAUC: A Python Package for Spatial Autocorrelation Analysis of Gene Sets

Spatial AUC is a powerful Python package designed to calculate spatial autocorrelation of gene sets using Moran's I statistic. It offers a streamlined and efficient way to analyze the spatial patterns of gene expression in spatial transcriptomics data.

## Features

- Retrieve gene sets from the MSigDB database
- Calculate area under the curve (AUC) for gene sets
- Compute spatial autocorrelation using Moran's I statistic
- Utilize Numba for performance optimization
- Seamlessly integrate with Scanpy and Squidpy for spatial analysis

## Installation

To install Spatial AUC, simply use pip:

```shell
pip install spatial-auc
```

Ensure that you have the following dependencies installed:

- numpy
- pandas
- scanpy
- squidpy
- numba
- gseapy

## Available Gene sets
For naming conventions and gene sets that can be retrieved consult the MSigDB website.
https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

Additionally, you can define your own gene sets and add them to the list.
The following structure is expected



## Usage

Here's a basic example of how to use Spatial AUC in your analysis:

```python
import scanpy as sc
from spatialAUC.spatialAUC import spatial_auc

# Load your Scanpy AnnData object
adata = sc.read_h5ad("your_data.h5ad")

# Calculate spatial autocorrelation
morans_table, adata_updated = spatial_auc(adata, gene_sets=["m5.all", "m2.all"], n_perms=1000)

# Access the results
print(morans_table)
```

For more detailed examples and usage instructions, please refer to the [Google Colab examples](https://colab.research.google.com/drive/1jWw7JBmPyL7-5L3U3ZvN4J4jZZZZ3Q3Q?usp=sharing).

## Large Dataset Optimization

If your dataset is too large for computations like this, consider using [Pseudovisium](https://github.com/BKover99/Pseudovisium) for hexagonal-binning. Pseudovisium can make your analysis an order of magnitude faster and more memory-efficient.

## Author

Bence Kover

- Twitter: [https://twitter.com/kover_bence](https://twitter.com/kover_bence)
- LinkedIn: [https://www.linkedin.com/in/ben-kover/](https://www.linkedin.com/in/ben-kover/)
