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
        highlight_thresh: float = None,
        species_taxa_file: str = None
) -> None:
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
   

    if species_taxa_file and highlight_thresh:
        with open(species_taxa_file, "r") as fh:
            species_taxa = [
                tuple(line.replace("\n", "").split("\t"))
                for line in fh.readlines()
            ]

        samples = zscores.columns.to_list()
        p_vals = highlight_df.iloc[:, 0]
        highlight_dict = {
            "x": list(), "y": list(),
            "tooltip": list(), "highlight": list(), "pair": list()
        }

        if os.isdir(highlight_data):
            files = os.listdir(highlight_data)
        else:
            files = highlight_data

        printf(f"Received files:\n{files}")
        for file in files:
            # TODO: will need to make this more general for other tables
            # 1) parameters?
            highlight_df = pd.read_csv(file, sep="\t").iloc[:, [1, 3, 8]]

            for i in range(len(p_vals)):
                if p_vals[i] < highlight_thresh:
                    sig_taxa = highlight_df.iloc[i, 0]
                    le_peps = highlight_df.iloc[i, 1].split("/")
                    for le_pep in le_peps:
                        highlight_dict["x"].append(
                            zscores.loc[le_pep, samples[0]]
                        )
                        highlight_dict["y"].append(
                            zscores.loc[le_pep, samples[1]]
                        )
                        highlight_dict["tooltip"].append(le_pep)
                        highlight_dict["highlight"].append(sig_taxa)
        for pair in pairs_sample_select:
            highlight_dict["pair"].append(f"{pair}")
        highlight_df = pd.DataFrame(highlight_dict)
        print(f"Highlight DataFrame:\n{highlight_df}")
        highlight_chart = alt.Chart(highlight_df).mark_point(
            filled=True, size=60
        ).encode(
            x=alt.X("x:Q", title=samples[0]),
            y=alt.Y("y:Q", title=samples[1]),
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
