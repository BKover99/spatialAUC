import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from numba import njit, prange
from gseapy import Msigdb
from tqdm import tqdm
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from typing import List, Union, Optional, Tuple



def get_df_from_gmt(
    categories: Union[List[str], str],
    version: str,
    genes: List[str],
    min_gene_ratio: float = 0.5,
    min_gene_count: int = 5,
    filter_go_kegg_reactome: bool = True
) -> pd.DataFrame:
    """
    Retrieve gene sets from MSigDB and return a DataFrame.

    This function fetches gene sets from the Molecular Signatures Database (MSigDB)
    based on specified categories and filters them according to various criteria.

    Args:
        categories: Categories of gene sets to retrieve. Can be a single category
            or a list of categories.
        version: Version of the MSigDB database to use.
        genes: List of genes specified by user or those present in the adata.var_names.
        min_gene_ratio: Minimum ratio of genes present in the gene set to the
            total genes in the gene set. Defaults to 0.5.
        min_gene_count: Minimum number of genes present in the gene set.
            Defaults to 5.
        filter_go_kegg_reactome: If True, filter gene sets to only include GO,
            KEGG, and REACTOME pathways. Defaults to True.

    Returns:
        A DataFrame containing filtered gene sets and their associated genes.

    Raises:
        ValueError: If no gene sets are found after filtering.
    """

    msig = Msigdb()
    if isinstance(categories, str):
        categories = [categories]
    dfs = []
    for category in categories:
        gmt = msig.get_gmt(category=category, dbver=version)
        df = pd.DataFrame(gmt.items(), columns=["gene_set", "genes"])
        dfs.append(df)
    if len(dfs) > 1:
        df = pd.concat(dfs, ignore_index=True)
    else:
        df = dfs[0]

    df["gene_present"] = df["genes"].apply(
        lambda x: len(set(x).intersection(set(genes)))
    )
    df["total_genes"] = df["genes"].apply(lambda x: len(x))
    df["gene_names_present"] = df["genes"].apply(
        lambda x: list(set(x).intersection(set(genes)))
    )
    df = df[
        (df["gene_present"] > min_gene_ratio * df["total_genes"])
        & (df["gene_present"] > min_gene_count)
    ]
    if filter_go_kegg_reactome:
        df = df[df["gene_set"].str.startswith(("GO", "KEGG", "REACTOME"))]
    df.reset_index(drop=True, inplace=True)
    return df


@njit(parallel=True)
def get_rank_numba(X: np.ndarray, axis: int = 0) -> np.ndarray:
    """
    Calculate the ranks of elements in a matrix along a specified axis using Numba.

    This function computes the ranks of values in the input matrix. It's optimized
    using Numba for parallel execution.

    Args:
        X: Input matrix.
        axis: Axis along which to calculate the ranks. 0 for rows, 1 for columns.
            Defaults to 0.

    Returns:
        A matrix of the same shape as the input, containing ranks.

    Note:
        This function is Numba-accelerated and designed for use within other functions.
    """
    if axis == 0:
        n_rows, n_cols = X.shape
        result = np.zeros_like(X, dtype=np.int64)
        for i in prange(n_rows):
            result[i] = np.argsort(np.argsort(-X[i]))
    else:
        n_rows, n_cols = X.shape
        result = np.zeros_like(X, dtype=np.int64)
        for i in prange(n_cols):
            result[:, i] = np.argsort(np.argsort(-X[:, i]))
    return result


def get_rank(adata: sc.AnnData, axis: int = 0) -> sc.AnnData:
    """
    Calculate the ranks of elements in an AnnData object along a specified axis.

    This function computes the ranks of gene expression values in the AnnData object.

    Args:
        adata: Input AnnData object.
        axis: Axis along which to calculate the ranks. 0 for rows (genes),
            1 for columns (cells). Defaults to 0.

    Returns:
        A new AnnData object with ranks calculated.

    Note:
        This function creates a copy of the input AnnData object to avoid
        modifying the original data.
    """

    adata = adata.copy()
    if not isinstance(adata.X, np.ndarray):
        adata.X = np.array(adata.X.todense())
    adata.X = get_rank_numba(adata.X, axis=axis)
    return adata


@njit
def calc_auc_numba(x_vals: np.ndarray, gene_indices: np.ndarray, axis: int = 0) -> float:
    """
    Calculate the area under the curve (AUC) using Numba.

    This function computes the AUC for a set of genes within a single cell or across cells.
    It's optimized using Numba for faster execution.

    Args:
        x_vals: Input values, typically ranks of genes.
        gene_indices: Indices of genes of interest within x_vals.
        axis: Axis along which to calculate the AUC. 0 for within a cell,
            1 for across cells. Defaults to 0.

    Returns:
        The calculated AUC value.

    Note:
        This function is Numba-accelerated and designed for use within other functions.
    """
    if axis == 0:
        y_vals = np.zeros(len(x_vals))
        indices = x_vals[gene_indices]
        for i in range(0, len(x_vals)):
            if i > 0:
                y_vals[i] = y_vals[i - 1]
            if i in indices:
                y_vals[i] += 1 / len(indices)
        return np.trapz(y_vals)
    else:
        y_vals = np.zeros(len(x_vals))
        x_vals = np.argsort(np.argsort(x_vals))
        indices = x_vals[gene_indices]
        for i in range(0, len(x_vals)):
            if i > 0:
                y_vals[i] = y_vals[i - 1]
            if i in indices:
                y_vals[i] += 1 / len(indices)
        return np.trapz(y_vals)


def get_auc(adata: sc.AnnData, df: pd.DataFrame, axis: int = 0) -> sc.AnnData:
    """
    Calculate the area under the curve (AUC) for gene sets in an AnnData object.

    This function computes AUC scores for each gene set in the provided DataFrame
    across all cells in the AnnData object.

    Args:
        adata: Input AnnData object.
        df: DataFrame containing gene sets and their associated genes.
        axis: Axis along which to calculate the AUC. 0 for rows (genes),
            1 for columns (cells). Defaults to 0.

    Returns:
        A new AnnData object with AUC values calculated. The `.X` attribute contains
        the AUC scores, with rows corresponding to cells and columns to gene sets.

    Note:
        This function uses Numba-accelerated functions for efficient calculation.
    
    """

    matrix = np.zeros((adata.shape[0], len(df)))
    adata = get_rank(adata, axis=axis)
    var_names = adata.var_names.tolist()

    for i in tqdm(range(len(df)), desc="Calculating AUC", unit="gene set"):
        genes = df.iloc[i]["genes"]
        gene_indices = np.array(
            [var_names.index(gene) for gene in genes if gene in var_names]
        )
        for j in range(adata.shape[0]):
            matrix[j, i] = calc_auc_numba(adata.X[j], gene_indices, axis=axis)

    new_adata = sc.AnnData(matrix)
    new_adata.obs = adata.obs
    new_adata.var = pd.DataFrame(df["gene_set"])
    new_adata.var_names = df["gene_set"]
    new_adata.obs_names = adata.obs_names
    if hasattr(adata, "obsp"):
        new_adata.obsp = adata.obsp
    if hasattr(adata, "obsm"):
        new_adata.obsm = adata.obsm
    if hasattr(adata, "uns"):
        new_adata.uns = adata.uns
    return new_adata



def compute_jaccard_similarities(contributing_genes: List[np.ndarray]) -> np.ndarray:
    """
    Compute pairwise Jaccard similarities between sets of genes.

    This function calculates the Jaccard similarity index between all pairs of gene sets
    provided in the input list.

    Args:
        contributing_genes: A list of numpy arrays, where each array contains gene names or indices.

    Returns:
        A square numpy array containing pairwise Jaccard similarities.

    Note:
        This function uses sparse matrices for efficient computation, making it suitable
        for large numbers of gene sets.
    """

    # Get unique genes across all arrays
    all_genes = np.unique(np.concatenate(contributing_genes))
    
    # Create a dictionary mapping each gene to a unique index
    gene_to_index = {gene: idx for idx, gene in enumerate(all_genes)}
    
    # Create a binary matrix where each row represents a set
    rows = []
    cols = []
    for i, genes in enumerate(contributing_genes):
        for gene in genes:
            rows.append(i)
            cols.append(gene_to_index[gene])
    
    # Create a sparse matrix
    data = np.ones(len(rows), dtype=int)
    matrix = csr_matrix((data, (rows, cols)), shape=(len(contributing_genes), len(all_genes)))
    
    # Compute Jaccard distances
    jaccard_distances = pdist(matrix.toarray(), metric='jaccard')
    
    # Convert distances to similarities
    jaccard_similarities = 1 - jaccard_distances
    
    # Convert to a square matrix
    similarity_matrix = squareform(jaccard_similarities)
    
    return similarity_matrix




def filter_gene_sets(new_morans_table: pd.DataFrame, contributing_genes: List[np.ndarray], similarity_threshold: float = 0.3) -> pd.DataFrame:
    """
    Filter gene sets based on Jaccard similarity.

    This function removes gene sets that are too similar to each other based on a
    Jaccard similarity threshold. It prioritizes keeping gene sets that appear earlier
    in the input DataFrame.

    Args:
        new_morans_table: DataFrame containing Moran's I statistics for gene sets.
        contributing_genes: List of arrays containing genes for each gene set.
        similarity_threshold: Jaccard similarity threshold above which gene sets
            are considered too similar. Defaults to 0.3.

    Returns:
        A filtered DataFrame with similar gene sets removed.

    Note:
        The order of gene sets in new_morans_table is assumed to represent their
        significance or priority.
    """

    
    # Compute Jaccard similarities
    similarities = compute_jaccard_similarities(contributing_genes)
    
    # Create a boolean mask for sets to keep
    keep_mask = np.ones(len(contributing_genes), dtype=bool)
    
    # Iterate through the sets in order of significance
    for i in range(len(contributing_genes)):
        if keep_mask[i]:
            # Find sets that are too similar to the current set
            similar_sets = np.where((similarities[i] > similarity_threshold) & (np.arange(len(contributing_genes)) > i))[0]
            
            # Mark these sets for removal
            keep_mask[similar_sets] = False
    
    # Apply the mask to the dataframe
    filtered_df = new_morans_table[keep_mask].copy()
    
    
    
    return filtered_df



def __spatial_auc__(
    adata: sc.AnnData,
    df: pd.DataFrame,
    n_perms: int = 1000,
    n_jobs: int = 2,
    axis: int = 0
) -> np.ndarray:
    """
    Internal function to calculate Moran's I for gene sets.

    This function computes AUC scores for gene sets and then calculates
    Moran's I spatial autocorrelation statistic for each gene set.

    Args:
        adata: Input AnnData object.
        df: DataFrame containing gene sets and their associated genes.
        n_perms: Number of permutations for calculating Moran's I p-values.
            Defaults to 1000.
        n_jobs: Number of jobs for parallel processing. Defaults to 2.
        axis: Axis along which to calculate AUC. 0 for rows (genes),
            1 for columns (cells). Defaults to 0.

    Returns:
        An array of Moran's I values for each gene set.

    Note:
        This is an internal function and should not be called directly by users.
    """

    adata = get_auc(adata, df, axis=axis)
    
    print("Calculating Moran's I for fake gene sets")
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=n_jobs)
    morans_I = adata.uns["moranI"]["I"].values

    return morans_I
    

import numpy as np
import pandas as pd

def generate_fake_gene_sets(contributing_genes: List[np.ndarray], n_fake_sets: int = 1000) -> pd.DataFrame:
    """
    Generate fake gene sets based on the characteristics of the original gene sets.

    This function creates random gene sets that mimic the size distribution and gene
    pool of the input gene sets. It's useful for creating null distributions for
    statistical testing.

    Args:
        contributing_genes: List of arrays containing genes for each original gene set.
        n_fake_sets: Number of fake gene sets to generate. Defaults to 1000.

    Returns:
        A DataFrame with columns 'gene_set' and 'genes', where 'genes' contains
        the randomly generated gene sets.

    Note:
        The generated fake gene sets maintain the size distribution of the original
        gene sets but with randomly selected genes.
    """
    
    # Get all unique genes from gene sets
    all_genes = np.unique(np.concatenate(contributing_genes))
    
    # Get lengths of original gene sets
    lengths = [len(x) for x in contributing_genes]
    
    # Draw n_perms lengths
    random_lengths = np.random.choice(lengths, n_fake_sets)
    
    # Draw random_lengths unique genes from all_genes
    random_sets = [np.random.choice(all_genes, length, replace=False) for length in random_lengths]
    
    # Create DataFrame with two columns "gene_set" and "genes"
    random_sets_df = pd.DataFrame({
        "gene_set": [f"random_set_{i}" for i in range(n_fake_sets)],
        "genes": random_sets
    })
    
    return random_sets_df


def spatial_auc(
    adata: sc.AnnData,
        df: Optional[pd.DataFrame] = None,
        genes: List[str] = [],
        gene_sets: List[str] = ["m5.all", "m2.all"],
        version: str = "2023.1.Mm",
        neighbors_defined: bool = True,
        n_perms: int = 1000,
        n_jobs: int = 2,
        axis: int = 0,
        min_gene_ratio: float = 0.3,
        min_gene_count: int = 5,
        sc: bool = False,
        jaccard_filter: float = 0.8,
        filter_go_kegg_reactome: bool = True,
        empirical_p: bool = True,
        n_fake_sets: int = 1000
    ) -> Tuple[pd.DataFrame, sc.AnnData]:
    """
    Calculate spatial autocorrelation of gene sets using Moran's I statistic.

    This function computes AUC scores for gene sets, calculates spatial
    autocorrelation using Moran's I, and optionally computes empirical p-values.

    Args:
        adata: Input AnnData object.
        df: Input DataFrame containing gene sets. If provided, it will be
            concatenated with the gene sets derived from MSigDB. Defaults to None.
        genes: List of genes to consider. If empty, all genes in adata.var_names
            will be used. Defaults to an empty list.
        gene_sets: List of gene set categories to retrieve from MSigDB.
            Defaults to ['m5.all', 'm2.all'].
        version: Version of the MSigDB database. Defaults to '2023.1.Mm'.
        neighbors_defined: Whether spatial neighbors have already been defined
            in the AnnData object. Defaults to True.
        n_perms: Number of permutations for calculating the p-value. Defaults to 1000.
        n_jobs: Number of jobs for parallel processing. Defaults to 2.
        axis: Rank genes within cells (0) or between cells (1). Defaults to 0.
        min_gene_ratio: Minimum ratio of genes present in the gene set to the
            total genes in the gene set. Defaults to 0.3.
        min_gene_count: Minimum number of genes present in the gene set. Defaults to 5.
        sc: If True, return the AnnData object without calculating spatial
            autocorrelation. Defaults to False.
        jaccard_filter: Jaccard similarity threshold for filtering gene sets.
            Defaults to 0.8.
        filter_go_kegg_reactome: If True, filter gene sets to only include GO,
            KEGG, and REACTOME pathways. Defaults to True.
        empirical_p: If True, calculate empirical p-values for the gene sets.
            Defaults to True.
        n_fake_sets: Number of fake gene sets to generate for empirical p-value
            calculation. Defaults to 1000.

    Returns:
        A tuple containing:
        - DataFrame with Moran's I statistics and p-values for each gene set.
        - Updated AnnData object with spatial gene set activity results.

    Raises:
        ValueError: If no gene sets are found after filtering.
    """

    print("Getting gene set")
    if genes == []:
        genes = adata.var_names.tolist()

    try:
        df_msigdb = get_df_from_gmt(
            gene_sets,
            version,
            genes,
            min_gene_ratio=min_gene_ratio,
            min_gene_count=min_gene_count,
            filter_go_kegg_reactome=filter_go_kegg_reactome,
        )
    except:
        print("Error: gene sets not found")
        df_msigdb = pd.DataFrame()

    if df is not None:
        try:
            if df.columns[0] == df_msigdb.columns[0]:
                df = pd.concat([df, df_msigdb], ignore_index=True)
            else:
                print("Error: df with wrong format")
                df = df_msigdb
        except:
            print("Error: df with wrong format")
            df = df_msigdb
    else:
        df = df_msigdb

    if df.empty:
        print("Error: no gene sets found")
        return None, adata

    print("Calculating AUC")
    adata_auc = get_auc(adata, df, axis=axis)

    if sc:
        return adata_auc
    else:
        if not neighbors_defined:
            print("Defining neighbors")
            sq.gr.spatial_neighbors(
                adata_auc, radius=250, coord_type="generic", delaunay=True
            )

        print("Calculating Moran's I")
        sq.gr.spatial_autocorr(adata_auc, mode="moran", n_perms=n_perms, n_jobs=n_jobs)
        morans_table_fr = adata_auc.uns["moranI"]

        contributing_genes = []
        for i in range(len(morans_table_fr)):
            contributing_genes.append(df_msigdb[df_msigdb["gene_set"] == morans_table_fr.index[i]].iloc[0]["gene_names_present"])
            
        morans_table_fr["contributing_genes"] = contributing_genes
        morans_table_fr["n_contributing_genes"] = [len(x) for x in contributing_genes]

        print("Filtering gene sets")
        #remove rows with nan
        morans_table_fr = morans_table_fr.dropna()
        morans_table_fr = filter_gene_sets(morans_table_fr, contributing_genes, similarity_threshold=jaccard_filter)
        
        
        #calculate empirical p_values for the gene sets
        if empirical_p:
            print("Generating fake gene sets")
            random_sets_df = generate_fake_gene_sets(contributing_genes, n_fake_sets=n_fake_sets)
            print("Calculating empirical p-values")
            empirical_morans_I = __spatial_auc__(adata, random_sets_df, n_perms=n_perms, n_jobs=n_jobs, axis=axis)
            morans_table_fr["p_empirical"] = [np.mean(empirical_morans_I > x) for x in morans_table_fr["I"]]


        return morans_table_fr, adata_auc
