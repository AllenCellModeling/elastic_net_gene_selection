# -*- coding: utf-8 -*-

import numpy as np
import altair as alt

from .genes import filter_out_unpenalized_genes
from .data import tidy

alt.data_transformers.enable("default", max_rows=None)


def thresh_lambda_plot(
    boot_results,
    adata,
    thresholds=np.linspace(0.01, 1, num=10),
    lambdas=np.geomspace(10, 0.01, num=10),
    unpenalized_genes=np.array([]),
):

    boot_betas = [
        filter_out_unpenalized_genes(
            beta=br["beta"],
            unpenalized_genes=unpenalized_genes,
            all_genes=[g.split("_")[0] for g in adata.var.index.values],
        )
        for br in boot_results
    ]
    gsel_bool = np.stack(boot_betas).mean(axis=0)
    gsel_thresh = [
        [(np.sum(gsel_bool[:, j] >= thresh)) for thresh in thresholds]
        for j, a in enumerate(lambdas)
    ]

    df_sel_thresh_alpha = tidy(np.stack(gsel_thresh))
    df_sel_thresh_alpha.columns = [
        "lambda index",
        "selection threshold index",
        "number of genes selected",
    ]
    df_sel_thresh_alpha["log number of genes selected"] = np.log1p(
        df_sel_thresh_alpha["number of genes selected"]
    )

    mychart = (
        alt.Chart(df_sel_thresh_alpha, width=600, height=600)
        .mark_rect()
        .encode(
            alt.X("selection threshold index:O", scale=alt.Scale(paddingInner=0)),
            alt.Y("lambda index:O", scale=alt.Scale(paddingInner=0)),
            alt.Color(
                "log number of genes selected:Q", scale=alt.Scale(scheme="greys")
            ),
        )
    )

    return mychart


def hub_persistence_plot(
    adata, boot_results, downsample_backgound=1000, width=400, height=400
):

    z = np.stack([br["beta"] for br in boot_results])
    sel_fracs = (z != 0).mean(axis=0)
    df_sel = tidy(sel_fracs)
    df_sel.columns = ["Gene Index", "Lambda Index", "Selection Fraction"]

    df_tmp = adata.var.copy()
    df_tmp["Gene Index"] = df_tmp.index.values.astype(np.int64)
    df_sel = df_sel.merge(df_tmp)
    df_sel = df_sel[df_sel["Gene Index"] < 1000 + downsample_backgound]

    chart = (
        alt.Chart(df_sel, width=width, height=height)
        .mark_line(opacity=0.25, interpolate="linear", color="blue")
        .encode(
            x="Lambda Index:Q",
            y="Selection Fraction:Q",
            detail="Gene Index:N",
            color="Type:N",
        )
    )

    return chart
