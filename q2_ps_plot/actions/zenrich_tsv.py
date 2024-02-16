#!/usr/bin/env python
from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSNPNFormat
)

import qiime2

# Name: zenrich_tsv
# Process: imports all TSV files into QZA files and runs the zenrich module
# with those inputs.
# Method inputs/parameters: ctx, data_filepath, zscores_filepath,
# negative_controls, negative_id, highlight_probes, source, pn_filepath,
# peptide_metadata, tooltip, color_by, negative_data_filepath, step_z_thresh,
# upper_z_thresh, lower_z_thresh, exact_z_thresh, exact_cs_thresh,
# pepsirf_binary
# Method outputs/Returned: the interactive visualization
# Dependencies: None
def zenrich_tsv(
        ctx,
        data_filepath,
        zscores_filepath,
        spline_x_filepath=None,
        spline_y_filepath=None,
        taxa_peps_md=None,  # TODO: change back to filepath
        p_val_thresh=0.05,
        p_vals=None,
        flex_reps=False,
        negative_controls=None,
        negative_id=None,
        source=None,
        pn_filepath=None,
        peptide_metadata=None,
        tooltip=["Species", "SpeciesID"],
        color_by="z_score_threshold",
        negative_data_filepath=None,
        highlighted_probes_filepath=None,
        step_z_thresh=5,
        upper_z_thresh=30,
        lower_z_thresh=5,
        exact_z_thresh=None,
        exact_cs_thresh="20",
        pepsirf_binary="pepsirf"
    ):

    # collect the zenrich module action
    zenrich = ctx.get_action("ps-plot", "zenrich")

    # import data into an artifact
    data = ctx.make_artifact(
        type="FeatureTable[Normed]",
        view=data_filepath,
        view_type=PepsirfContingencyTSVFormat
    )

    # import zscores into an artifact
    zscores = ctx.make_artifact(
        type="FeatureTable[Zscore]",
        view=zscores_filepath,
        view_type=PepsirfContingencyTSVFormat
    )
    
    # if negative data provided import into an artifact
    if negative_data_filepath:
        negative_data = ctx.make_artifact(
            type="FeatureTable[Normed]",
            view=negative_data_filepath,
            view_type=PepsirfContingencyTSVFormat
        )
    # otherwisde set negative data as none
    else:
        negative_data = None

    # if highlihgted probes provided import into an aftifact
    if highlighted_probes_filepath:
        highlighted_probes = ctx.make_artifact(
            type="InfoSNPN",
            view=highlighted_probes_filepath,
            view_type=PepsirfInfoSNPNFormat
        )
    # otherwise set highlighed probes to none
    else:
        highlighted_probes = None

    # run the zenrich module with all the inputs and parameters given
    zenrich_vis, = zenrich(
        data=data,
        zscores=zscores,
        spline_x_filepath=spline_x_filepath,
        spline_y_filepath=spline_y_filepath,
        taxa_peps_md=taxa_peps_md,
        p_val_thresh=p_val_thresh,
        p_vals=p_vals,
        flex_reps=flex_reps,
        negative_controls=negative_controls,
        negative_id=negative_id,
        highlight_probes=highlighted_probes,
        source=source,
        pn_filepath=pn_filepath,
        peptide_metadata=peptide_metadata,
        tooltip=tooltip,
        color_by=color_by,
        negative_data=negative_data,
        step_z_thresh=step_z_thresh,
        upper_z_thresh=upper_z_thresh,
        lower_z_thresh=lower_z_thresh,
        exact_z_thresh=exact_z_thresh,
        exact_cs_thresh=exact_cs_thresh,
        pepsirf_binary=pepsirf_binary
    )

    #return the visualization
    return zenrich_vis
