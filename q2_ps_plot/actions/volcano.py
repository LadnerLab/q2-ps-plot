#!/usr/bin/env python
import altair as alt
import os
import numpy as np
import pandas as pd


# TODO: suggest making this action as general as possible
def volcano(
        output_dir: str, 
        x: list = None,
        y: list = None,
        taxa: list = None,
        xy_dir: str = None,
        xy_access: list = ["x", "y"],
        taxa_access: str = None,
        x_threshold: float = 0.4,
        y_thresholds: list = [0.05],
        log: bool = True,
        xy_labels: list = ["x", "y"],
        titles: list = [""]
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.disable_max_rows()

    if (xy_dir is not None
            and x is None
            and y is None
            and taxa is None
            and xy_access is not None
            and taxa_access is not None):
        x = []
        y = []
        taxa = []
        files = os.listdir(xy_dir)
        files = sorted(files)
        for file in files:
            if "~" not in file:
                continue
            data = pd.read_csv(f"{xy_dir}/{file}", sep="\t")
            x.append(data.loc[:, xy_access[0]].to_list())
            y.append(data.loc[:, xy_access[1]].to_list())
            taxa.append(data.loc[:, taxa_access].to_list())
        titles = sorted(titles)
    elif (xy_dir is None
            and x is not None
            and y is not None
            and taxa is not None):
        x = [x]
        y = [y]
        taxa = [taxa]
    else:
        print(
            "A directory with tables, as well as unrelated parameters were"
            "provided:\n"
            f"x = {x}\n"
            f"y = {y}\n"
            f"taxa = {taxa}\n"
            f"xy_access = {xy_access}\n"
            f"taxa_access = {taxa_access}\n"
            "Please review and edit parameters accordingly."
        )
        return

    charts = []
    for i in range(len(x)):
        titles_len = len(titles)
        if log:
            log_adjusted_y = -np.log10(y[i])
            volcano_dict = { "y": log_adjusted_y, "x": x[i] }
            sig_taxa_df = make_sig_taxa_df(
                log_adjusted_y, x[i], taxa[i],
                y[i], y_thresholds[i], x_threshold
            )
            sort = "ascending"
        else:
            volcano_dict = {"y": y[i], "x": x[i]}
            sig_taxa_df = make_sig_taxa_df(
                y[i], x[i], taxa[i],
                y[i], y_thresholds[i], x_threshold
            )
            sort = "descending"
        volcano_df = pd.DataFrame(volcano_dict)

        volcano_chart = alt.Chart(
            volcano_df, title=titles[i % titles_len]
        ).mark_circle(
            size=50, color="black"
        ).encode(
            x=alt.X(
                "x:Q",
                title=xy_labels[0]
            ),
            y=alt.Y(
                "y:Q",
                title=xy_labels[1],
                sort=sort
            )
        )
        chart = alt.layer(volcano_chart)

        if sig_taxa_df is not None:
            sig_taxa = alt.Chart(
                sig_taxa_df, title=titles[i % titles_len]
            ).mark_circle(size=60).encode(
                x=alt.X(
                    "x:Q",
                    title=xy_labels[0]
                ),
                y=alt.Y(
                    "y:Q",
                    title=xy_labels[1],
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
            chart = alt.layer(chart, sig_taxa).resolve_scale(
                color="independent"
            )
        charts.append(chart)

    final_chart = alt.vconcat(*charts)
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
