import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from numba import njit, prange
from gseapy import Msigdb
from tqdm import tqdm


def get_df_from_gmt(
    categories,
    version,
    genes,
    min_gene_ratio=0.5,
    min_gene_count=5,
    filter_go_kegg_reactome=True,
):
    """
    Retrieve gene sets from MSigDB and return a DataFrame.

    Args:
        categories (list or str): Categories of gene sets to retrieve.
        version (str): Version of the MSigDB database.
        genes (list): List of genes present in the adata.var_names.
        min_gene_ratio (float): Minimum ratio of genes present in the gene set to the total genes in the gene set.
        min_gene_count (int): Minimum number of genes present in the gene set.
        filter_go_kegg_reactome (bool): If True, filter gene sets to only include GO, KEGG, and REACTOME pathways.

    Returns:
        pandas.DataFrame: DataFrame containing gene sets and their associated genes.
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
def get_rank_numba(X, axis=0):
    """
    Calculate the ranks of elements in a matrix along a specified axis using Numba.

    Args:
        X (numpy.ndarray): Input matrix.
        axis (int): Axis along which to calculate the ranks (0 for rows, 1 for columns).

    Returns:
        numpy.ndarray: Matrix of ranks.
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


def get_rank(adata, axis=0):
    """
    Calculate the ranks of elements in an AnnData object along a specified axis.

    Args:
        adata (scanpy.AnnData): Input AnnData object.
        axis (int): Axis along which to calculate the ranks (0 for rows, 1 for columns).

    Returns:
        scanpy.AnnData: AnnData object with ranks calculated.
    """
    adata = adata.copy()
    if not isinstance(adata.X, np.ndarray):
        adata.X = np.array(adata.X.todense())
    adata.X = get_rank_numba(adata.X, axis=axis)
    return adata


@njit
def calc_auc_numba(x_vals, gene_indices, axis=0):
    """
    Calculate the area under the curve (AUC) using Numba.

    Args:
        x_vals (numpy.ndarray): Input values.
        gene_indices (numpy.ndarray): Indices of genes.
        axis (int): Axis along which to calculate the AUC (0 for rows, 1 for columns).

    Returns:
        float: AUC value.
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


def get_auc(adata, df, axis=0):
    """
    Calculate the area under the curve (AUC) for gene sets in an AnnData object.

    Args:
        adata (scanpy.AnnData): Input AnnData object.
        df (pandas.DataFrame): DataFrame containing gene sets and their associated genes.
        axis (int): Axis along which to calculate the AUC (0 for rows, 1 for columns).

    Returns:
        scanpy.AnnData: AnnData object with AUC values calculated.
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


def spatial_auc(
    adata,
    df=None,
    genes=[],
    gene_sets=["m5.all", "m2.all"],
    version="2023.1.Mm",
    neighbors_defined=True,
    n_perms=1000,
    n_jobs=2,
    axis=0,
    min_gene_ratio=0.3,
    min_gene_count=5,
    sc=False,
    filter_go_kegg_reactome=True,
):
    """
    Calculate spatial autocorrelation of gene sets using Moran's I statistic.

    Args:
        adata (scanpy.AnnData): Input AnnData object.
        df (pandas.DataFrame, optional): Input DataFrame containing gene sets. If provided, it will be concatenated with the gene sets derived from MSigDB. Defaults to None.
        genes (list): List of genes to consider. If empty, all genes in adata.var_names will be used. Defaults to [].
        gene_sets (list): List of gene set categories to retrieve from MSigDB. Defaults to ['m5.all', 'm2.all'].
        version (str): Version of the MSigDB database. Defaults to '2023.1.Mm'.
        neighbors_defined (bool): Whether spatial neighbors have already been defined in the AnnData object. Defaults to True.
        n_perms (int): Number of permutations for calculating the p-value. Defaults to 1000.
        n_jobs (int): Number of jobs for parallel processing. Defaults to 2.
        axis (int): Rank genes within cells (0) or between cells (1). Defaults to 0, which is the normal AUC calculation.
        min_gene_ratio (float): Minimum ratio of genes present in the gene set to the total genes in the gene set. Defaults to 0.3.
        min_gene_count (int): Minimum number of genes present in the gene set. Defaults to 5.
        sc (bool): If True, return the AnnData object without calculating spatial autocorrelation. Defaults to False.
        filter_go_kegg_reactome (bool): If True, filter gene sets to only include GO, KEGG, and REACTOME pathways. Defaults to True.

    Returns:
        tuple: A tuple containing two elements:
            - morans_table_fr (pandas.DataFrame): DataFrame containing Moran's I statistics and p-values for each gene set.
            - adata (scanpy.AnnData): Updated AnnData object with spatial autocorrelation results.
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
    adata = get_auc(adata, df, axis=axis)

    if sc:
        return adata
    else:
        if not neighbors_defined:
            print("Defining neighbors")
            sq.gr.spatial_neighbors(
                adata, radius=250, coord_type="generic", delaunay=True
            )

        print("Calculating Moran's I")
        sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=n_jobs)
        morans_table_fr = adata.uns["moranI"]

        return morans_table_fr, adata
