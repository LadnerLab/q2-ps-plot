import altair as alt
import numpy as np
import os
import pandas as pd
import qiime2
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import math

from q2_pepsirf.format_types import PepsirfContingencyTSVFormat


def zscatter(
        output_dir: str,
        zscores: PepsirfContingencyTSVFormat,
        pairs_file: str,
        colors_file: str = "",
        p_val_access: str = None,
        le_peps_access: str = None,
        taxa_access: str = None,
        spline_file: str = None,
        highlight_data: str = None,
        highlight_threshold: float = 0.05,
        vis_outputs_dir: str = None
) -> None:
    alt.data_transformers.disable_max_rows()

    if vis_outputs_dir:
        plot_output_dir = os.path.join(vis_outputs_dir, "scatter_plots")
        os.mkdir(plot_output_dir)

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

    unsorted_pairs = list()
    pair_2_title = dict()
    with open(pairs_file, "r") as fh:
        fh.readline()
        for line in fh.readlines():
            line_tup = tuple(line.replace("\n", "").split("\t"))
            pair = line_tup[0:2]
            unsorted_pairs.append(f"{pair[0]}~{pair[1]}")

            if( len(line_tup) > 2):
                pair_2_title[f"{pair[0]}~{pair[1]}"] = line_tup[2]
            else:
                pair_2_title[f"{pair[0]}~{pair[1]}"] = ""

    pairs = sorted(unsorted_pairs)

    if highlight_data and not os.path.isfile(highlight_data):
        files = os.listdir(highlight_data)
        path = highlight_data
        files = sorted(files)
    else:
        files = highlight_data
        path = "."

    sample_dropdown = alt.binding_select(options=unsorted_pairs, name="Sample Select")
    sample_select = alt.selection_point(
        fields=["pair"],
        bind=sample_dropdown,
        name="pair",
        value=[{"pair": unsorted_pairs[0]}]
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
            "peptide": list(), "taxa": list(), "pair": list()
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
                        highlight_dict["peptide"].append(le_pep)
                        highlight_dict["taxa"].append(sig_taxa)
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
        spline_df = pd.DataFrame(spline_df)
        spline_chart = alt.Chart(
            spline_df
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

    # create color scale
    color_scale=alt.Scale(range=[
                    "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00",
                    "#CC79A7"])
    shape=alt.Shape("taxa:N", legend=None)
    legend=alt.Legend(title="Significant Taxa")

    if colors_file:
        all_species = list(set(highlight_df["taxa"].to_list()))
        num_extra_colors = len(all_species)

        color_df = pd.read_csv(colors_file, sep="\t", header=None, names=["Taxa", "Color"])
        # remove species that are not in all_species
        color_df = color_df[color_df["Taxa"].isin(all_species)]
        species_list = color_df["Taxa"].to_list()
        colors_list = color_df["Color"].to_list()
        num_extra_colors -= len(species_list)

        color_iter = iter(plt.cm.rainbow(np.linspace(0, 1, num_extra_colors)))

        # fill lists for species not included in color file
        for species in all_species:
            if species not in species_list:
                species_list.append(species)
                colors_list.append(clr.to_hex(next(color_iter)))

        color_scale=alt.Scale(domain=species_list, range=colors_list)
        shape=alt.Shape("taxa:N", scale=alt.Scale(domain=species_list), legend=None)
        legend=alt.Legend(title="Significant Taxa", columns=int(math.ceil(len(species_list)/30)), symbolLimit=0)

    if highlight_df is not None:
        highlight_chart = alt.Chart(highlight_df).mark_point(
            filled=True, size=60
        ).encode(
            x=alt.X("x:Q"),
            y=alt.Y("y:Q"),
            color=alt.Color(
                "taxa:N",
                scale=color_scale,
                legend=legend
            ),
            # https://github.com/altair-viz/altair/issues/1181
            shape=shape,
            tooltip=["peptide", "taxa"]
        ).transform_filter(
            sample_select
        )
        final_chart = alt.layer(final_chart, highlight_chart).resolve_scale(
            color="independent",
            shape="independent"
        )

    titleDf = pd.DataFrame(pair_2_title.items(), columns=["pair", "title"])
    title = alt.Chart(titleDf).mark_text(
        size=30, dx=chart_width/2
    ).encode(
        text='title:N'
    ).transform_filter(
        sample_select
    )
    final_chart = alt.vconcat(title, final_chart)
    
    final_chart.save(os.path.join(output_dir, "index.html"))

    if vis_outputs_dir:
        heatmap_pair_df = {pair: df for pair, df in heatmap_df.groupby("pair")}
        if spline_file:
            spline_pair_df = {pair: df for pair, df in spline_df.groupby("pair")}
        if highlight_data:
            highlight_pair_df = {pair: df for pair, df in highlight_df.groupby("pair")}

        
        for pair in pairs:
            final_chart = None
            if pair in heatmap_pair_df.keys():
                heatmap_chart = alt.Chart(
                    heatmap_pair_df[pair], width=chart_width, height=chart_height, 
                    title=alt.TitleParams(pair, anchor='middle')
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
                )
                final_chart = heatmap_chart
                
            if spline_file and pair in spline_pair_df.keys():
                spline_chart = alt.Chart(
                    spline_pair_df[pair]
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
                if final_chart != None:
                    final_chart = alt.layer(final_chart, spline_chart)
                else:
                    final_chart = spline_chart

            if highlight_df is not None and pair in highlight_pair_df.keys():
                highlight_chart = alt.Chart(highlight_pair_df[pair]).mark_point(
                    filled=True, size=60
                ).encode(
                    x=alt.X("x:Q"),
                    y=alt.Y("y:Q"),
                    color=alt.Color(
                        "taxa:N",
                        scale=color_scale,
                        legend=legend
                    ),
                    # https://github.com/altair-viz/altair/issues/1181
                    shape=shape,
                    tooltip=["peptide", "taxa"]
                )
                if final_chart != None:
                    final_chart = alt.layer(final_chart, highlight_chart).resolve_scale(
                        color="independent",
                        shape="independent"
                    )
                else:
                    final_chart = highlight_chart
            
            if final_chart != None:
                final_chart.save(os.path.join(plot_output_dir, f"{pair}_scatter.html"))
            else:
                print(f"Skipped scatter plot for {pair}")