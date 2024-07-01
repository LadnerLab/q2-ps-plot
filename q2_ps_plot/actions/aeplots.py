#!/usr/bin/env python

import altair as alt
import pandas as pd
import os

def aeplots(
        output_dir: str, 
        pos_nes_ae_file: str,
        neg_nes_ae_file: str,
        xy_access: list = ["Events", "Species"],
        xy_labels: list = ["Number of AEs in cohort", "Species"]
) -> None:
    alt.data_transformers.disable_max_rows()

    pos_df = pd.read_csv(pos_nes_ae_file, sep="\t")
    neg_df = pd.read_csv(neg_nes_ae_file, sep="\t")

    pos_df["NES"]="Positive"
    neg_df["NES"]="Negative"

    ae_df = pd.concat([pos_df, neg_df])

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
            f"{xy_access[1]}:N"
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

    final_chart.save(os.path.join(output_dir, "index.html"));