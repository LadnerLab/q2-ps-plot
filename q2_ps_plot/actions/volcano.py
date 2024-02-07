#!/usr/bin/env python
import altair as alt
import os
import pandas as pd


def volcano(
    output_dir: str, 
    p_vals: list,
    es: list,
    p_val_thresh: float,
    es_thresh: float
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.enable("default", max_rows=None)

    volcano_dict = {
        col_name: list() for col_name in ["p-vals", "es"]
    }
    volcano_dict["p-vals"] = p_vals
    volcano_dict["es"] = es
    volcano_df = pd.DataFrame(volcano_dict)

    volcano_chart = alt.Chart(volcano_df).mark_circle(size=50).encode(
        x=alt.X(
            "es:Q",
            title="Enrichment score"
        ),
        y=alt.Y(
            "p-vals:Q",
            title="Adjusted p-values",
            sort="descending"
        ),
        # TODO: integrate tooltip
    )

    sig_taxa_dict = {["id"]: list(), ["p-vals"]: list(), ["es"]}
    sig_taxa_df = pd.DataFrame(sig_taxa_dict)
    
    sig_taxa = alt.Chart(sig_taxa_df).mark_circle(size=60).encode(
        x=alt.X(
            "es:Q",
            title="Enrichment score"
        ),
        y=alt.Y(
            "p-vals:Q",
            title="p-values",
            sort="descending"
        ),
        # TODO: integrate color
    )

    final_chart = alt.layer(
        volcano_chart,
        sig_taxa
    ).properties(
        title="PSEA"  # TODO: get suggestion for a better name
    ).resolve_scale(color="independent")

    final_chart.save(os.path.join(output_dir, "index.html"))
