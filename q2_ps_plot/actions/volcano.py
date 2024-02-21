#!/usr/bin/env python
import altair as alt
import os
import numpy as np
import pandas as pd


# TODO: suggest making this action as general as possible
def volcano(
        output_dir: str, 
        y: list,
        x: list,
        taxa: list,
        y_thresh: float = 0.05,
        x_thresh: float = 0.4,
        log: bool = True,
        x_label: str = "x",
        y_label: str = "y"
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.enable("default", max_rows=None)

    if log:
        log_adjusted_y = -np.log10(y)
        volcano_dict = {"y": log_adjusted_y, "x": x}
        sig_taxa_df = make_sig_taxa_df(
            log_adjusted_y, x, taxa, y, y_thresh, x_thresh
        )
        sort = "ascending"
    else:
        volcano_dict = {"y": y, "x": x}
        sig_taxa_df = make_sig_taxa_df(y, x, taxa, y, y_thresh, x_thresh)
        sort = "descending"
    volcano_df = pd.DataFrame(volcano_dict)

    volcano_chart = alt.Chart(volcano_df).mark_circle(
        size=50, color="black"
    ).encode(
        x=alt.X(
            "x:Q",
            title=x_label
        ),
        y=alt.Y(
            "y:Q",
            title=y_label,
            sort=sort
        ),
    )
    final_chart = alt.layer(volcano_chart)

    if sig_taxa_df is not None:
        sig_taxa = alt.Chart(sig_taxa_df).mark_circle(size=60).encode(
            x=alt.X(
                "x:Q",
                title=x_label
            ),
            y=alt.Y(
                "y:Q",
                title=y_label,
                sort=sort
            ),
            color=alt.Color(
                "taxa:N",
                scale=alt.Scale(range=[
                    "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00",
                    "#CC79A7"
                ]),
                legend=None
            ),
            tooltip="taxa"
        )
        final_chart = alt.layer(final_chart, sig_taxa).resolve_scale(
            color="independent"
        )

    final_chart.save(os.path.join(output_dir, "index.html"))


def make_sig_taxa_df(
    y,
    x,
    taxa,
    y_test,
    y_thresh,
    x_thresh
) -> pd.DataFrame:
    if taxa is None:
        return None

    sig_taxa_dict = {"taxa": list(), "y": list(), "x": list()}
    for i in range(len(taxa)):
        if y_test[i] < y_thresh and abs(x[i]) > x_thresh:
            sig_taxa_dict["taxa"].append(taxa[i])
            sig_taxa_dict["y"].append(y[i])
            sig_taxa_dict["x"].append(x[i])
    return pd.DataFrame(sig_taxa_dict)
