#!/usr/bin/env python
import altair as alt
import os
import pandas as pd


def volcano(
    output_dir: str, 
    p_vals: list,
    es: list
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.enable("default", max_rows=None)

    volcano_dict = {
        col_name: list() for col_name in ["p_vals", "es"]
    }
    volcano_dict["p_vals"] = p_vals
    volcano_dict["es"] = es
    volcano_df = pd.DataFrame(volcano_dict)

    volcano_chart = alt.Chart(volcano_df).mark_circle(size=50).encode(
        x=alt.X(
            "es:Q",
            title="Enrichment score"
        ),
        y=alt.Y(
            "p_vals:Q",
            title="Adjusted p-values"
        ),
        # TODO: integrate tooltip
    )

    final_chart = alt.layer(volcano_chart).properties(title="PSEA")  # TODO: get suggestion for a better name

    final_chart.save(os.path.join(output_dir, "index.html"))
