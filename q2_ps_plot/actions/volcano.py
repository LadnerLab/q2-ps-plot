#!/usr/bin/env python
import altair as alt
import os
import pandas as pd


# TODO: suggest making this action as general as possible
def volcano(
    output_dir: str, 
    # TODO: this could be a metadata table?
    p_vals: list,
    es: list,
    taxa: list,
    p_val_thresh: float,
    es_thresh: float = 0.4
) -> None:
    # TODO: change into temp directory
    alt.data_transformers.enable("default", max_rows=None)

    volcano_dict = {"p-vals": p_vals, "es": es}
    volcano_df = pd.DataFrame(volcano_dict)

    volcano_chart = alt.Chart(volcano_df).mark_circle(
        size=50, color="black"
    ).encode(
        x=alt.X(
            "es:Q",
            title="Enrichment score"
        ),
        y=alt.Y(
            "p-vals:Q",
            title="p-values",
            sort="descending"
        ),
        # TODO: integrate tooltip?
    )

    # TODO: if this is as general as possible, should this be asked?
    if taxa:
        sig_taxa_dict = {"id": list(), "p-vals": list(), "es": list()}
        # TODO: assuming taxa, p_vals, and es are the same length
        for i in range(len(taxa)):
            if p_vals[i] < p_val_thresh and abs(es[i]) > es_thresh:
                sig_taxa_dict["id"].append(taxa[i])
                sig_taxa_dict["p-vals"].append(p_vals[i])
                sig_taxa_dict["es"].append(es[i])
        sig_taxa_df = pd.DataFrame(sig_taxa_dict)
        print(f"Significant Taxa DataFrame:\n{sig_taxa_df}")
    
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
    else:
        # TODO: create final chart without highlighting
        pass

    final_chart.save(os.path.join(output_dir, "index.html"))
