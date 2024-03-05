import altair as alt
import numpy as np
import os
import pandas as pd
import qiime2

from q2_pepsirf.format_types import PepsirfContingencyTSVFormat


def zscatter(
        output_dir: str,
        zscores: PepsirfContingencyTSVFormat,
        pairs_file: str,
        spline_file: str = None,
        highlight_data: str = None,
        highlight_thresholds: list = None,
        species_taxa_file: str = None
) -> None:
    alt.data_transformers.disable_max_rows()

    zscores = zscores.view(pd.DataFrame)
    zscores = zscores.transpose()

    spline_data_df = None
    if spline_file:
        spline_data_df = pd.read_csv(spline_file, sep="\t")

    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    pairs = sorted(pairs)
    with open(species_taxa_file, "r") as fh:
        species_taxa = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]

    if not os.path.isfile(highlight_data):
        files = os.listdir(highlight_data)
        path = highlight_data
    else:
        files = highlight_data
        path = "."
    files = sorted(files)

    charts = []
    f = 0
    for file in files:
        if "~" not in file:
            continue
        pair = pairs[f]
        print(f"File: {file}\nPair: {pair}\n")

        heatmap_dict = {
            "bin_x_start": list(), "bin_x_end": list(),
            "bin_y_start": list(), "bin_y_end": list(), "count": list()
        }
    
        x = zscores.loc[:, pair[0]]
        y = zscores.loc[:, pair[1]]
        heatmap, x_edges, y_edges = np.histogram2d(x, y, bins=(70, 70))
        for x in range(0, heatmap.shape[0]):
            for y in range(0, heatmap.shape[1]):
                # assumes count is for a peptide
                count = heatmap[x, y]
                if count == 0.0:
                    continue
                bin_x_start = x_edges[x]
                bin_x_end = x_edges[x + 1]
                bin_y_start = y_edges[y]
                bin_y_end = y_edges[y + 1]

                heatmap_dict["bin_x_start"].append(bin_x_start)
                heatmap_dict["bin_x_end"].append(bin_x_end)
                heatmap_dict["bin_y_start"].append(bin_y_start)
                heatmap_dict["bin_y_end"].append(bin_y_end)
                heatmap_dict["count"].append(count)
        heatmap_df = pd.DataFrame(heatmap_dict)
        xy_max = heatmap_df.loc[:, ["bin_x_end", "bin_y_end"]].max()
        ratio = (xy_max[0] / xy_max[1]) - 1
        chart_height = 500
        chart_width = chart_height + (20*ratio)

        heatmap_chart = alt.Chart(
            heatmap_df, width=chart_width, height=chart_height
        ).mark_rect().encode(
            alt.X("bin_x_start:Q", title=pair[0]),
            alt.X2("bin_x_end:Q", title=None),
            alt.Y("bin_y_start:Q", title=pair[1]),
            alt.Y2("bin_y_end:Q", title=None),
            alt.Color(
                "count:Q",
                scale=alt.Scale(scheme="greys"),
                legend=alt.Legend(title="Point Frequency")
            )
        )
        chart = alt.layer(heatmap_chart)
        
        
        if spline_data_df is not None:
            spline_dict = {
                "x": spline_data_df.loc[:, pair[0]],
                "y": spline_data_df.loc[:, pair[1]]
            }
        
            spline_df = pd.DataFrame(spline_dict)
            spline_chart = alt.Chart(
                pd.DataFrame(spline_df)
            ).mark_square(size=20).encode(
                x=alt.X("x:Q"),
                y=alt.Y("y:Q"),
                color=alt.Color(
                    "x:N",
                    scale=alt.Scale(range=["#FF0000"]),
                    # reference: https://github.com/altair-viz/altair/issues/620
                    legend=None
                )
            )
            chart = alt.layer(chart, spline_chart)   


        if species_taxa_file and highlight_thresholds:
            samples = zscores.columns.to_list()
            highlight_dict = {
                "x": list(), "y": list(),
                "tooltip": list(), "highlight": list()
            }

            print(f"Working with highlight file: {path}/{file}")
            highlight_df = pd.read_csv(
                f"{path}/{file}", sep="\t"
            ).iloc[:, [3, 4, 8]]
            p_vals = highlight_df.iloc[:, 0].to_list()

            for i in range(len(p_vals)):
                if p_vals[i] < highlight_thresholds[f]:
                    sig_taxa = highlight_df.iloc[i, 2]
                    le_peps = highlight_df.iloc[i, 1].split("/")
                    for le_pep in le_peps:
                        highlight_dict["x"] = zscores.loc[le_pep, samples[0]]
                        highlight_dict["y"] = zscores.loc[le_pep, samples[1]]
                        highlight_dict["tooltip"].append(le_pep)
                        highlight_dict["highlight"].append(sig_taxa)
            highlight_df = pd.DataFrame(highlight_dict)
            highlight_chart = alt.Chart(highlight_df).mark_point(
                filled=True, size=60
            ).encode(
                x=alt.X("x:Q", title=pair[0]),
                y=alt.Y("y:Q", title=pair[1]),
                color=alt.Color(
                    "highlight:N",
                    scale=alt.Scale(range=[
                        "#E69F00", "#56B4E9", "#009E73",
                        "#F0E442", "#0072B2", "#D55E00",
                        "#CC79A7"
                    ]),
                    legend=alt.Legend(title="Significant Taxa")
                ),
                # https://github.com/altair-viz/altair/issues/1181
                shape=alt.Shape(
                    "highlight:N",
                    legend=None
                ),
                tooltip="tooltip"
            )
            chart = alt.layer(chart, highlight_chart).resolve_scale(
                color="independent",
                shape="independent"
            )

        charts.append(chart)
        f += 1

    final_chart = alt.concat(*charts)
    final_chart.save(os.path.join(output_dir, "index.html"))
    fnal
