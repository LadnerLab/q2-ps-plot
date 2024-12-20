#!/usr/bin/env python

import pandas as pd
import altair as alt
import os
import qiime2

def compareCS_histogram(
    output_dir: str,
    sample_names: list,
    category_names: list,
    c_counts: list,
    positions: list,
    num_bins: int
):
    chart_dict = {
        "Sample Name": sample_names,
        "Category Name": category_names,
        "C count": c_counts,
        "Position": positions
    }

    chart_df = pd.DataFrame(chart_dict)

    # create sample dropdown menu
    samples = list(chart_df["Sample Name"].unique())
    sample_dropdown = alt.binding_select(
        options=samples, name="Sample Select"
    )
    sample_select = alt.selection_point(
        fields=["Sample Name"], bind=sample_dropdown,
        name="Sample Name", value=[{"Sample Name": samples[0]}]
    )

    category_groups = chart_df.groupby("Category Name")

    charts = list()
    for category, category_chart_df in category_groups:
        charts.append(
            alt.Chart(category_chart_df).mark_bar(opacity=0.75).encode(
                alt.X("Position:Q").bin(maxbins=num_bins),
                alt.Y("C count:Q"),
                alt.Color("Category Name:N")
            ).transform_filter(
                sample_select
            )
        )


    final_chart = alt.layer(
            *charts
        ).add_params(
            sample_select
        )

    final_chart.save(os.path.join(output_dir, "index.html"), scale_factor=10.0)