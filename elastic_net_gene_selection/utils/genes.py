# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import swifter  # noqa F401


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


def get_thresh_lambda_df(
    boot_results,
    gene_names=[],
    thresholds=np.linspace(0.01, 1, num=100),
    lambdas=np.geomspace(10, 0.01, num=100),
):
    """
    get names of selected genes for each combination of thresholds and lambdas.
    returns a df of results, with the list of gene names merged into one string 
    for eachthreshold/lambda combo.

    TODO: make this work (efficiently) with unpenalized genes
    """

    # fraction of times each gene is selected at each lambda
    gsel_fracs = np.stack([b["beta"] for b in boot_results]).mean(axis=0)

    # make sure we're all lined up
    assert gsel_fracs.shape == (len(gene_names), len(lambdas))

    # make a df: one row per threshold + lambda (index) combo
    index = pd.MultiIndex.from_product(
        [thresholds, range(len(lambdas))], names=["selection threshold", "lambda index"]
    )
    df_thresh_lam = pd.DataFrame(index=index).reset_index()

    # we iterate on the index of each lambda, but really we want the lambdas in there
    df_thresh_lam["lambda"] = df_thresh_lam["lambda index"].map(
        dict(enumerate(lambdas))
    )

    # apply a name-getting-function to each row in parallel (using swifter)
    df_thresh_lam["selected genes"] = df_thresh_lam.swifter.apply(
        lambda row: gene_names[
            gsel_fracs[:, int(row["lambda index"])] > row["selection threshold"]
        ],
        axis="columns",
    )

    df_thresh_lam["number of genes selected"] = df_thresh_lam.swifter.apply(
        lambda row: len(row["selected genes"]), axis="columns",
    )

    df_thresh_lam["selected genes"] = df_thresh_lam["selected genes"].str.join(", ")

    # drop the temporary lambda index column
    df_thresh_lam = df_thresh_lam.drop(["lambda index"], axis="columns")

    return df_thresh_lam
