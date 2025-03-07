#!/usr/bin/env python

import altair as alt
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np

def aeplots(
        output_dir: str, 
        pos_nes_ae_file: str,
        neg_nes_ae_file: str,
        colors_file: str = "",
        xy_access: list = ["Events", "Species"],
        xy_labels: list = ["Number of AEs in cohort", "Species"],
        vis_outputs_dir: str = None
) -> None:
    alt.data_transformers.disable_max_rows()

    pos_df = pd.read_csv(pos_nes_ae_file, sep="\t")
    neg_df = pd.read_csv(neg_nes_ae_file, sep="\t")

    pos_df["NES"]="Positive"
    neg_df["NES"]="Negative"

    ae_df = pd.concat([pos_df, neg_df])

    # create color scale
    color_scale=alt.Scale(range=[
                    "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00",
                    "#CC79A7"])
    if colors_file:
        all_species = list(set(ae_df["Species"].to_list()))
        num_extra_colors = len(all_species)

        color_df = pd.read_csv(colors_file, sep="\t", header=None, names=["Species", "Color"])
        # remove species that are not in all_species
        color_df = color_df[color_df["Species"].isin(all_species)]
        species_list = color_df["Species"].to_list()
        colors_list = color_df["Color"].to_list()
        num_extra_colors -= len(species_list)

        color_iter = iter(plt.cm.rainbow(np.linspace(0, 1, num_extra_colors)))

        # fill lists for species not included in color file
        for species in all_species:
            if species not in species_list:
                species_list.append(species)
                colors_list.append(clr.to_hex(next(color_iter)))

        color_scale=alt.Scale(domain=species_list, range=colors_list)

    bar_chart = alt.Chart(ae_df).mark_bar().encode(
        alt.X(
            f"{xy_access[0]}:Q", 
            title=xy_labels[0]
            ),
        alt.Y(
            f"{xy_access[1]}:N",
            axis=alt.Axis(grid=True),
            sort=alt.EncodingSortField(field=xy_access[1], op="count", order='ascending'),
            title=None
            ),
        alt.Color(
            f"{xy_access[1]}:N",
            scale=color_scale,
            legend=None
            )
        )

    text = bar_chart.mark_text(
        align='left',
        baseline='middle',
        dx=3
        ).encode(
            text=f"{xy_access[0]}:Q"
        )

    final_chart = (bar_chart + text).facet(
            column=alt.Column(
                "NES:N", 
                sort=alt.SortField(field='NES', order='descending'), 
                title=None,
                header=alt.Header(labelFontSize=15, labelFontWeight='bold')
                )
        ).resolve_scale(
            y='independent'
        )

    final_chart.save(os.path.join(output_dir, "index.html"))
    if vis_outputs_dir:
        final_chart.save(os.path.join(vis_outputs_dir, "aeplots.html"))