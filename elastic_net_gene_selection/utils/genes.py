# -*- coding: utf-8 -*-

import numpy as np


def get_gene_set(df, num_genes=25):
    gene_subset_sizes = df["Number of genes"].unique()
    closest_gene_subset_size = gene_subset_sizes[
        np.abs(gene_subset_sizes - num_genes).argmin()
    ]
    closest_gene_subset = df[df["Number of genes"] == closest_gene_subset_size][
        "Gene name"
    ].values
    return closest_gene_subset


def filter_out_unpenalized_genes(beta, unpenalized_genes, all_genes):
    beta_out = beta.copy()
    assert beta_out.shape[0] == len(all_genes)
    unpenalized_genes_bool_inds = np.array(
        [gene in set(unpenalized_genes) for gene in all_genes]
    )
    beta_out[unpenalized_genes_bool_inds, :] = False
    return beta_out


def get_selected_genes(
    boot_results,
    adata,
    lambda_index=75,
    selection_threshold_index=75,
    thresholds=np.linspace(0.01, 1, num=10),
    unpenalized_genes=np.array([]),
):
    if len(unpenalized_genes) > 0:
        boot_betas = [
            filter_out_unpenalized_genes(
                beta=br["beta"],
                unpenalized_genes=unpenalized_genes,
                all_genes=adata.var.index,
            )
            for br in boot_results
        ]
    else:
        boot_betas = [br["beta"] for br in boot_results]

    gsel_bool = np.stack(boot_betas).mean(axis=0)
    genes_bool = gsel_bool[:, lambda_index] > thresholds[selection_threshold_index]
    return adata.var["gene_name"][genes_bool].values.astype(str)
