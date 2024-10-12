# spatialAUC: A Python Package for Spatial Autocorrelation Analysis of Gene Sets

spatialAUC is a powerful Python package designed to calculate spatial autocorrelation of gene sets using Moran's I statistic. It offers a streamlined and efficient way to analyze the spatial patterns of gene expression in spatial transcriptomics data. AUC might stand for autocorrelation, or area under the curve, or both.

<img width="381" alt="Screenshot 2024-06-21 at 19 15 03" src="https://github.com/BKover99/spatialAUC/assets/91386576/5784116d-958c-4188-95c2-e2edfd9865cb">
<img width="586" alt="Screenshot 2024-06-21 at 19 14 59" src="https://github.com/BKover99/spatialAUC/assets/91386576/d63cbb5d-4aea-44b8-a444-6632a1ea9060">


## Features

- Retrieve gene sets from the MSigDB database
- Calculate area under the curve (AUC) for gene sets
- Compute spatial autocorrelation using Moran's I statistic
- Utilize Numba for performance optimization
- Seamlessly integrate with Scanpy and Squidpy for spatial analysis

## Installation

To install spatialAUC, simply use pip:

```shell
pip install spatialauc
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
The following structure is expected:
 
 ### Columns
 
 - **gene_set**: The name of the gene set, formatted as a string.
 - **genes**: A list of all genes associated with the gene set, formatted as a list of strings.

 
 ### Example Table
 
 | gene_set                       | genes                                                                 |
 |--------------------------------|-----------------------------------------------------------------------|
 | GOBP_NEUROTRANSMITTER_UPTAKE   | [Slc6a5, Slc18a1, Slc6a3, Slc29a2, Drd2, Drd3, ...]                   |
 | GOBP_FEVER_GENERATION          | [Ccr5, Cnr1, Ednrb, Il1a, Il1b, Il1rn, Ptger3, ...]                   |



## Usage

Here's a basic example of how to use spatialAUC in your analysis:

```python
import scanpy as sc
from spatialAUC.spatialAUC import spatial_auc

# Load your Scanpy AnnData object
adata = sc.read_h5ad("your_data.h5ad")

#Processing
#... normalise, call spatial_neighbors, etc...

# Calculate spatial autocorrelation - all args displayed here
morans_table, adata_updated = spatial_auc(adata, df=None, genes=[], gene_sets=['m5.all', 'm2.all'],
                                         version='2023.1.Mm', neighbors_defined=True, n_perms=1000,
                                         n_jobs=2,axis=0, min_gene_ratio=0.3, min_gene_count=5, sc=False,
                                         jaccard_filter = 0.8, filter_go_kegg_reactome=True,
                                         empirical_p = True, n_fake_sets = 1000)

# Access the results
print(morans_table)
```

For more detailed examples and usage instructions, please refer to the examples below:
[SpatialAUC tutorial on Xenium whole mouse pup data](https://github.com/BKover99/spatialAUC/blob/main/Tutorials/spatialAUC_tutorial_Xenium_whole_mouse_pup.ipynb)
[SpatialAUC tutorial on CosMx mouse brain data](https://github.com/BKover99/spatialAUC/blob/main/Tutorials/spatialAUC_tutorial_cosmx_mouse_brain.ipynb)

## Empirical P-values

SpatialAUC calculates empirical p-values to enhance the statistical robustness of its spatial autocorrelation analysis. This process involves:

- Generating 'fake' gene sets that mimic the size distribution and gene pool of the original sets.
- Calculating Moran's I for these fake sets to create a null distribution.
- Determining the empirical p-value for each real gene set as the proportion of fake sets with Moran's I values exceeding that of the real set.

This data-driven approach accounts for dataset-specific characteristics and reduces potential biases in traditional p-value calculations. By default, spatialAUC uses 1000 fake gene sets, but users can adjust this number to balance computational cost with statistical precision. Empirical p-values provide a more reliable measure of statistical significance, especially in complex spatial transcriptomics datasets.

## Removing Redundant Gene Sets

SpatialAUC implements a filtering mechanism to remove redundant gene sets, ensuring that the final results are not skewed by highly similar sets. This process uses the Jaccard similarity index to measure the overlap between gene sets. The algorithm iterates through the gene sets in order of their Moran's I statistics. For each set, it calculates its Jaccard similarity with all subsequent sets. If the similarity exceeds a user-defined threshold (default is 0.8), the less significant set is removed. This approach prioritizes keeping gene sets with stronger spatial autocorrelation while eliminating redundancies. Users can adjust the similarity threshold to control the stringency of this filtering process. By removing redundant gene sets, spatialAUC provides a more concise and interpretable set of results, focusing on distinct biological pathways and processes that show significant spatial patterns.


## Large Dataset Optimization - Pooling cells by hexagonal-binning

If your dataset is too large for computations like this (or just want to pool nearby cells), consider using [Pseudovisium](https://github.com/BKover99/Pseudovisium) for hexagonal-binning. Pseudovisium can make your analysis an order of magnitude faster and more memory-efficient.

## Citation
To cite spatialAUC, cite the Pseudovisium pre-print.
https://www.biorxiv.org/content/10.1101/2024.07.23.604776v1

DOI: https://doi.org/10.1101/2024.07.23.604776

*Kover B, Vigilante A. Rapid and memory-efficient analysis and quality control of large spatial transcriptomics datasets [Internet]. bioRxiv; 2024. p. 2024.07.23.604776. Available from: https://www.biorxiv.org/content/10.1101/2024.07.23.604776v1*


## Author

Bence Kover

- Twitter: [https://twitter.com/kover_bence](https://twitter.com/kover_bence)
- LinkedIn: [https://www.linkedin.com/in/ben-kover/](https://www.linkedin.com/in/ben-kover/)
