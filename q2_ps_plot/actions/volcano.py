#!/usr/bin/env python

import altair as alt
import os
import numpy as np
import pandas as pd


def volcano(
        output_dir: str, 
        x: list = None,
        y: list = None,
        taxa: list = [""],
        xy_dir: str = None,
        xy_access: list = ["x", "y"],
        taxa_access: str = None,
        x_threshold: float = 0.4,
        y_threshold: float = 0.05,
        log: bool = True,
        xy_labels: list = ["x", "y"],
        titles: list = [""]
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.disable_max_rows()

    if xy_dir:
        # TODO: make assertions here?
        x = []
        y = []
        files = os.listdir(xy_dir)
        files = sorted(files)
        for file in files:
            if "~" not in file:
                continue
            data = pd.read_csv(f"{xy_dir}/{file}", sep="\t")
            x.append(data.loc[:, xy_access[0]].to_list())
            y.append(data.loc[:, xy_access[1]].to_list())
            taxa.append(data.loc[:, taxa_access].to_list())
    elif x and y:
        x = [x]
        y = [y]
        taxa = [taxa]
    titles = sorted(titles)

    volcano_dict = { "x": list(), "y": list(), "pair": list() }

    for i in range(len(x)):
        if log:
            log_adjusted_y = -np.log10(y[i])
            volcano_dict["y"].extend(log_adjusted_y)
            # TODO: highlight significant taxa
            sort = "ascending"
        else:
            volcano_dict["y"].extend(y[i])
            # TODO: highlight significant taxa
            sort = "descending"
        volcano_dict["x"].extend(x[i])
        volcano_dict["pair"].extend([titles[i]] * len(x[i]))
    volcano_df = pd.DataFrame(volcano_dict)
    volcano_df.to_csv("volcano_df.tsv", sep="\t")

    sample_dropdown = alt.binding_select(options=titles, name="Sample Select")
    sample_select = alt.selection_point(
        fields=["pair"],
        bind=sample_dropdown,
        name="pair",
        value=[{"pair": titles[0]}]
    )

    volcano_chart = alt.Chart(volcano_df).mark_circle(
        size=50, color="black"
    ).encode(
        x=alt.X("x:Q", title=xy_labels[0]),
        y=alt.Y("y:Q", title=xy_labels[1], sort=sort)
    ).add_params(
        sample_select
    ).transform_filter(
        sample_select
    )
    final_chart = alt.layer(volcano_chart)

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

    sig_taxa_dict = { "y": list(), "x": list(), "taxa": list() }
    for i in range(len(taxa)):
        if y_test[i] < y_thresh and abs(x[i]) > x_thresh:
            sig_taxa_dict["y"].append(y[i])
            sig_taxa_dict["x"].append(x[i])
            sig_taxa_dict["taxa"].append(taxa[i])
    return pd.DataFrame(sig_taxa_dict)
