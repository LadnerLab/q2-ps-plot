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
        p_val_access: str = None,
        le_peps_access: str = None,
        taxa_access: str = None,
        spline_file: str = None,
        highlight_data: str = None,
        highlight_threshold: float = 0.05
) -> None:
    alt.data_transformers.disable_max_rows()

    if highlight_data:
        assert p_val_access, \
            f"'highlight_data' was provided, but nothing was provided for" \
            " 'p_val_access'!"
        assert le_peps_access, \
            f"'highlight_data' was provided, but nothing was provided for" \
            " 'le_peps_access'!"
        assert taxa_access, \
            f"'highlight_data' was provided, but nothing was provided for" \
            " 'taxa_access'!"

    zscores = zscores.view(pd.DataFrame)
    zscores = zscores.transpose()

    with open(pairs_file, "r") as fh:
        pairs = [
            line.replace("\n", "").replace("\t", "~")
            for line in fh.readlines()
        ]
    pairs = sorted(pairs)

    if highlight_data and not os.path.isfile(highlight_data):
        files = os.listdir(highlight_data)
        path = highlight_data
        files = sorted(files)
    else:
        files = highlight_data
        path = "."

    sample_dropdown = alt.binding_select(options=pairs, name="Sample Select")
    sample_select = alt.selection_point(
        fields=["pair"],
        bind=sample_dropdown,
        name="pair",
        value=[{"pair": pairs[0]}]
    )

    heatmap_dict = {
        "bin_x_start": list(), "bin_x_end": list(),
        "bin_y_start": list(), "bin_y_end": list(),
        "count": list(), "pair": list()
    }
    p = 0
    for pair in pairs:
        pair = pair.split("~")

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
                heatmap_dict["pair"].append(pairs[p])
        p += 1
    heatmap_df = pd.DataFrame(heatmap_dict)
    xy_max = heatmap_df.loc[:, ["bin_x_end", "bin_y_end"]].max()
    ratio = (xy_max[0] / xy_max[1]) + 1
    chart_height = 500
    chart_width = chart_height + (20*ratio)

    highlight_df = None
    if highlight_data:
        highlight_dict = {
            "x": list(), "y": list(),
            "tooltip": list(), "highlight": list(), "pair": list()
        }

        f = 0
        for file in files:
            pair = pairs[f].split("~")
            if "~" not in file and pairs[f] not in file.rsplit(".", 1)[0]:
                continue
            highlight_df = pd.read_csv(
                f"{path}/{file}", sep="\t"
            ).loc[:, [p_val_access, le_peps_access, taxa_access]]

            p_vals = highlight_df.iloc[:, 0].to_list()

            for i in range(len(p_vals)):
                if p_vals[i] < highlight_threshold:
                    le_peps = highlight_df.iloc[i, 1].split("/")
                    sig_taxa = highlight_df.iloc[i, 2]
                    for le_pep in le_peps:
                        highlight_dict["x"].append(
                            zscores.loc[le_pep, pair[0]]
                        )
                        highlight_dict["y"].append(
                            zscores.loc[le_pep, pair[1]]
                        )
                        highlight_dict["tooltip"].append(le_pep)
                        highlight_dict["highlight"].append(sig_taxa)
                        highlight_dict["pair"].append(pairs[f])
            f += 1
        highlight_df = pd.DataFrame(highlight_dict)

    heatmap_chart = alt.Chart(
        heatmap_df, width=chart_width, height=chart_height
    ).mark_rect().encode(
        alt.X("bin_x_start:Q", title="Time Point 1"),
        alt.X2("bin_x_end:Q"),
        alt.Y("bin_y_start:Q", title="Time Point 2"),
        alt.Y2("bin_y_end:Q"),
        alt.Color(
            "count:Q",
            scale=alt.Scale(scheme="greys"),
            legend=alt.Legend(title="Point Frequency")
        )
    ).add_params(
        sample_select
    ).transform_filter(
        sample_select
    )
    final_chart = alt.layer(heatmap_chart)
        
    if spline_file:
        spline_df = pd.read_csv(spline_file, sep="\t")
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
        ).transform_filter(
            sample_select
        )
        final_chart = alt.layer(final_chart, spline_chart)

    if highlight_df is not None:
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
            shape=alt.Shape(
                "highlight:N",
                legend=None
            ),
            tooltip="tooltip"
        ).transform_filter(
            sample_select
        )
        final_chart = alt.layer(final_chart, highlight_chart).resolve_scale(
            color="independent",
            shape="independent"
        )

    final_chart.save(os.path.join(output_dir, "index.html"))
