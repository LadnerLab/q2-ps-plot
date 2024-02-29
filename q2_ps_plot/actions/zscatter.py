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
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    alt.data_transformers.enable("default", max_rows=None)
    zscores = zscores.view(pd.DataFrame)
    zscores = zscores.transpose()

    heatmap_dict = {
        "bin_x_start": list(), "bin_x_end": list(), "bin_y_start": list(),
        "bin_y_end": list(), "count": list(), "pair": list()
    }

    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    
    pairs_sample_select = []
    for pair in pairs:
        x = zscores.loc[:, pair[0]]
        y = zscores.loc[:, pair[1]]
        pair_sample_select = f"{pair[0]}~{pair[1]}"
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
                heatmap_dict["pair"].append(pair_sample_select)
        pairs_sample_select.append(pair_sample_select)
    heatmap_df = pd.DataFrame(heatmap_dict)
    xy_max = heatmap_df.loc[:, ["bin_x_end", "bin_y_end"]].max()
    ratio = (xy_max[0] / xy_max[1]) - 1
    chart_height = 500
    chart_width = chart_height + (20*ratio)

    pairs_binding = alt.binding_select(
        options=pairs_sample_select, name="Pair select"
    )
    pairs_select = alt.selection_point(
        fields=["pair"],
        bind=pairs_binding,
        name="pair",
        value=[{"pair": pairs_sample_select[0]}]
    )

    heatmap_chart = alt.Chart(
        heatmap_df, width=chart_width, height=chart_height
    ).mark_rect().encode(
        alt.X("bin_x_start:Q"),
        alt.X2("bin_x_end:Q"),
        alt.Y("bin_y_start:Q"),
        alt.Y2("bin_y_end:Q"),
        alt.Color(
            "count:Q",
            scale=alt.Scale(scheme="greys"),
            legend=alt.Legend(title="Point Frequency")
        )
    ).add_params(pairs_select).transform_filter(pairs_select)
    final_chart = alt.layer(heatmap_chart)
   

    if species_taxa_file and highlight_thresholds:
        samples = zscores.columns.to_list()
        highlight_dict = {
            "x": list(), "y": list(),
            "tooltip": list(), "highlight": list(), "pair": list()
        }

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
            path = ""

        print(f"Received files:\n{files}")
        h_thresh_idx = 0
        for file in files:
            if "~" not in file:
                continue
            # TODO: will need to make this more general for other tables
            # 1) parameters?
            highlight_df = pd.read_csv(
                f"{path}/{file}", sep="\t"
            ).iloc[:, [3, 4, 8]]
            p_vals = highlight_df.iloc[:, 0].to_list()

            print(f"len(p_vals) = {len(p_vals)}\nlen(highlight_thresholds) = {len(highlight_thresholds)}\n")
            x = []
            y = []
            le_peps = []
            sig_taxas = []
            for i in range(len(p_vals)):
                if p_vals[i] < highlight_thresholds[h_thresh_idx]:
                    sig_taxas.append(highlight_df.iloc[i, 2])
                    curr_le_peps = highlight_df.iloc[i, 1].split("/")
                    print(f"curr_le_peps: {curr_le_peps}")
                    print(f"x = {zscores.loc[curr_le_peps, samples[0]]}")
                    print(f"y = {zscores.loc[curr_le_peps, samples[1]]}")
                    x.append(zscores.loc[curr_le_peps, samples[0]].to_list())
                    y.append(zscores.loc[curr_le_peps, samples[1]].to_list())
                    le_peps.append(curr_le_peps)
            highlight_dict["x"].append(x)
            highlight_dict["y"].append(y)
            highlight_dict["tooltip"].append(le_peps)
            highlight_dict["highlight"].append(sig_taxas)
            h_thresh_idx += 1

        for pair in pairs_sample_select:
            highlight_dict["pair"].append(f"{pair}")
        for k, v in highlight_dict.items():
            print(f"{k}:")
            for i in v:
                print(f"\t(type: {type(i)}) = {i}")
        highlight_df = pd.DataFrame(highlight_dict)
        highlight_chart = alt.Chart(highlight_df).mark_point(
            filled=True, size=60
        ).encode(
            x=alt.X("x:Q"),
            y=alt.Y("y:Q"),
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
            # shape=alt.Shape(
            #     "highlight:N",
            #     legend=None
            # ),
            tooltip="tooltip"
        ).transform_filter(pairs_select)
        final_chart = alt.layer(
            final_chart,
            highlight_chart
        ).resolve_scale(
            color="independent",
            shape="independent"
        )


    # if spline_file:
    #     spline_df = pd.read_csv(spline_file, sep="\t")
    #     spline_chart = alt.Chart(
    #         pd.DataFrame(spline_df)
    #     ).mark_square(size=20).encode(
    #         x="x:Q",
    #         y="y:Q",
    #         color=alt.Color(
    #             "x:N",
    #             scale=alt.Scale(range=["#FF0000"]),
    #             # reference: https://github.com/altair-viz/altair/issues/620
    #             legend=None
    #         )
    #     )
    #     final_chart = alt.layer(final_chart, spline_chart)

    final_chart.save(os.path.join(output_dir, "index.html"))


def fill_highlight_dict():
    pass
