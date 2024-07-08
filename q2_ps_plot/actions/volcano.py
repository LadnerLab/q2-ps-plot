#!/usr/bin/env python

import altair as alt
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import math

def volcano(
        output_dir: str, 
        pairs_file: str,
        colors_file: str = "",
        x: list = None,
        y: list = None,
        taxa: list = None,
        xy_dir: str = None,
        xy_access: list = ["x", "y"],
        taxa_access: str = None,
        x_threshold: float = 0.4,
        y_threshold: float = 0.05,
        log: bool = True,
        xy_labels: list = ["x", "y"]
) -> None:
    alt.data_transformers.disable_max_rows()

    if xy_dir:
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
            if taxa_access:
                taxa.append(data.loc[:, taxa_access].to_list())
    elif x and y:
        x = [x]
        y = [y]
        if not taxa:
            taxa = []
        else:
            taxa = [taxa]
    unsorted_pairs = list()
    pair_2_title = dict()
    with open(pairs_file, "r") as fh:
        fh.readline()
        for line in fh.readlines():
            line_tup = tuple(line.replace("\n", "").split("\t"))
            pair = line_tup[0:2]
            unsorted_pairs.append(f"{pair[0]}~{pair[1]}")

            if( len(line_tup) > 2):
                pair_2_title[f"{pair[0]}~{pair[1]}"] = line_tup[2]
            else:
                pair_2_title[f"{pair[0]}~{pair[1]}"] = ""

    pairs = sorted(unsorted_pairs)

    sample_dropdown = alt.binding_select(options=unsorted_pairs, name="Sample Select")
    sample_select = alt.selection_point(
        fields=["pair"],
        bind=sample_dropdown,
        name="pair",
        value=[{"pair": unsorted_pairs[0]}]
    )

    volcano_dict = { "x": list(), "y": list(), "pair": list() }

    charting_ys = []
    for i in range(len(x)):
        if log:
            charting_ys.append(-np.log10(y[i]))
            sort = "ascending"
        else:
            charting_ys.append(y[i])
            sort = "descending"
        volcano_dict["x"].extend(x[i])
        volcano_dict["y"].extend(charting_ys[i])
        volcano_dict["pair"].extend([pairs[i]] * len(x[i]))
    volcano_df = pd.DataFrame(volcano_dict)

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

    highlight_df = None
    if taxa:
        highlight_dict = {
            "x": list(), "y": list(),
            "taxa": list(), "pair": list()
        }
        for i in range(len(taxa)):
            for j in range(len(taxa[i])):
                if y[i][j] < y_threshold and abs(x[i][j]) > x_threshold:
                    highlight_dict["x"].append(x[i][j])
                    highlight_dict["y"].append(charting_ys[i][j])
                    highlight_dict["taxa"].append(taxa[i][j])
                    highlight_dict["pair"].append(pairs[i])
        highlight_df = pd.DataFrame(highlight_dict)

        # create color scale
        color_scale=alt.Scale(range=[
                    "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00",
                    "#CC79A7"])
        legend=alt.Legend(title="Significant Taxa")
        if colors_file:
            all_species = list(set(highlight_df["taxa"].to_list()))
            num_extra_colors = len(all_species)

            color_df = pd.read_csv(colors_file, sep="\t", header=None, names=["Taxa", "Color"])
            # remove species that are not in all_species
            color_df = color_df[color_df["Taxa"].isin(all_species)]
            species_list = color_df["Taxa"].to_list()
            colors_list = color_df["Color"].to_list()
            num_extra_colors -= len(species_list)

            color_iter = iter(plt.cm.rainbow(np.linspace(0, 1, num_extra_colors)))

            # fill lists for species not included in color file
            for species in all_species:
                if species not in species_list:
                    species_list.append(species)
                    colors_list.append(clr.to_hex(next(color_iter)))

            color_scale=alt.Scale(domain=species_list, range=colors_list)
            legend=alt.Legend(title="Significant Taxa", columns=int(math.ceil(len(species_list)/30)), symbolLimit=0)

        highlight_chart = alt.Chart(highlight_df).mark_circle(
            size=60, filled=True, opacity=1.0
        ).encode(
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
                scale=color_scale,
                legend=legend
            ),
            tooltip="taxa"
        ).add_params(
            sample_select
        ).transform_filter(
            sample_select
        )
        final_chart = alt.layer(final_chart, highlight_chart).resolve_scale(
            color="independent"
        )

    titleDf = pd.DataFrame(pair_2_title.items(), columns=["pair", "title"])
    title = alt.Chart(titleDf).mark_text(
        size=25, dx=150
    ).encode(
        text='title:N'
    ).transform_filter(
        sample_select
    )
    final_chart = alt.vconcat(title, final_chart)

    final_chart.save(os.path.join(output_dir, "index.html"))
