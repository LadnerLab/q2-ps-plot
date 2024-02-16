#!/usr/bin/env python
import altair as alt
import os
import numpy as np
import pandas as pd


# TODO: suggest making this action as general as possible
def volcano(
        output_dir: str, 
        p_vals: list,
        es: list,
        taxa: list,
        x_label: str = None,
        y_label: str = None,
        p_val_thresh: float = 0.05,
        es_thresh: float = 0.4
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.enable("default", max_rows=None)

    log_adjusted_p_vals = -np.log10(p_vals)

    volcano_dict = {"adjust-p-vals": log_adjusted_p_vals, "es": es}
    volcano_df = pd.DataFrame(volcano_dict)

    volcano_chart = alt.Chart(volcano_df).mark_circle(
        size=50, color="black"
    ).encode(
        x=alt.X(
            "es:Q",
            title=x_label
        ),
        y=alt.Y(
            "adjust-p-vals:Q",
            title=y_label
        ),
        # TODO: integrate tooltip?
    )

    sig_taxa_dict = {"id": list(), "adjust-p-vals": list(), "es": list()}
    # TODO: assuming taxa, log_adjusted_p_vals, and es are the same length
    for i in range(len(taxa)):
        if p_vals[i] < p_val_thresh and abs(es[i]) > es_thresh:
            sig_taxa_dict["id"].append(taxa[i])
            sig_taxa_dict["adjust-p-vals"].append(log_adjusted_p_vals[i])
            sig_taxa_dict["es"].append(es[i])
    sig_taxa_df = pd.DataFrame(sig_taxa_dict)

    sig_taxa = alt.Chart(sig_taxa_df).mark_circle(size=60).encode(
        x=alt.X(
            "es:Q",
            title=x_label
        ),
        y=alt.Y(
            "adjust-p-vals:Q",
            title=y_label
        ),
        color=alt.Color(
            "id:N",
            scale=alt.Scale(range=[
                "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00",
                "#CC79A7"
            ]),
            legend=None  # TODO: do we want a legend?
        ),
        tooltip="id"
    )

    final_chart = alt.layer(
        volcano_chart,
        sig_taxa
    ).properties(
        title="PSEA"  # TODO: get suggestion for a better name
    ).resolve_scale(color="independent")

    final_chart.save(os.path.join(output_dir, "index.html"))
