#!/usr/bin/env python

import pandas as pd
import altair as alt
import os

def compareCS_scatter(
    output_dir: str,
    x: list,
    y: list,
    c_codenames: list,
    s_codenames: list,
    parent_sequences: list,
    c_counts: list,
    sample_names: list,
    color_scheme_1: str = "dark2",
    color_scheme_2: str = "blues"
):
    # make sure all lists are the same length
    assert len(x)==len(y) and len(y)==len(c_codenames) and len(c_codenames)==len(s_codenames) and len(s_codenames)==len(sample_names), \
        "Input lists are not all the same length!"
    
    chart_dict = {
        "C z-score": x,
        "S z-score": y,
        "C codename": c_codenames,
        "S codename": s_codenames,
        "Parent Sequence": parent_sequences,
        "C count": c_counts,
        "Sample Name": sample_names
    }

    max_value = max(x+y)

    chart_df = pd.DataFrame(chart_dict)

    chart_1 = alt.Chart(chart_df).mark_circle(size=50).encode(
        alt.X(
            "C z-score:Q",
            title="Cysteine version z-score",
            scale=alt.Scale(domain=[0,max_value], nice=False)
            ),
        alt.Y(
            "S z-score:Q",
            title="Serine version z-score",
            scale=alt.Scale(domain=[0,max_value], nice=False)
            ),
        alt.Color(
            "Sample Name:N",
            scale=alt.Scale(
                scheme=color_scheme_1
                )
            ),
        tooltip=[
            "C codename:N",
            "S codename:N",
            "C z-score:Q",
            "S z-score:Q",
            "Parent Sequence:N",
            "Sample Name:N"
            ]
        )
    
    chart_2 = alt.Chart(chart_df).mark_circle(size=50).encode(
        alt.X(
            "C z-score:Q",
            title="Cysteine version z-score",
            scale=alt.Scale(domain=[0,max_value], nice=False)
            ),
        alt.Y(
            "S z-score:Q",
            title="Serine version z-score",
            scale=alt.Scale(domain=[0,max_value], nice=False)
            ),
        alt.Color(
            "C count:Q",
            scale=alt.Scale(
                scheme=color_scheme_2
                )
            ),
        tooltip=[
            "C codename:N",
            "S codename:N",
            "C z-score:Q",
            "S z-score:Q",
            "Parent Sequence:N",
            "Sample Name:N"
            ]
        )
    
    line = pd.DataFrame({
    'x': [0, max_value],
    'y': [0, max_value],
    })

    line_plot = alt.Chart(line).mark_line(
        color='gray',
        strokeDash=[5, 5],
        strokeWidth=1
    ).encode(
        x='x',
        y='y'
    )

    scatter_plot_colored_by_species = line_plot + chart_1
    scatter_plot_colored_by_c_count = line_plot + chart_2

    final_plot = alt.vconcat(scatter_plot_colored_by_species, scatter_plot_colored_by_c_count)

    final_plot = final_plot.configure_view(
                            continuousHeight=500,
                            continuousWidth=500
                            )

    final_plot.save(os.path.join(output_dir, "index.html"))
