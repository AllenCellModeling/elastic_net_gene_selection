# -*- coding: utf-8 -*-

import numpy as np
import altair as alt

from .data import tidy

alt.data_transformers.enable("default", max_rows=None)


def thresh_lambda_plot(df):
    mychart = (
        alt.Chart(df, width=600, height=600)
        .mark_rect()
        .encode(
            x=alt.X(
                "selection threshold:O",
                scale=alt.Scale(paddingInner=0),
                axis=alt.Axis(format=".2f"),
            ),
            y=alt.Y(
                "lambda:O", scale=alt.Scale(paddingInner=0), axis=alt.Axis(format=".2f")
            ),
            color=alt.Color(
                "number of genes selected:Q",
                scale=alt.Scale(scheme="greys", type="symlog"),
            ),
            tooltip=["selected genes"],
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
