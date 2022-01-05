#!/usr/bin/env python
import pandas as pd
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat, PepsirfInfoSNPNFormat
import qiime2

def zenrich_tsv(ctx,
            data_filepath,
            zscores_filepath,
            negative_controls,
            source = None,
            pn_filepath = None,
            peptide_metadata = None,
            tooltip = ['Species', 'SpeciesID'],
            color_by = 'z_score_threshold',
            negative_data_filepath = None,
            highlighted_probes_filepath = None,
            step_z_thresh=5,
            upper_z_thresh=30,
            lower_z_thresh=5,
            exact_z_thresh=None,
            pepsirf_binary="pepsirf"):

        zenrich = ctx.get_action('ps-plot', 'zenrich')

        data = ctx.make_artifact(type='FeatureTable[Normed]',
                                        view=data_filepath,
                                        view_type=PepsirfContingencyTSVFormat)

        zscores = ctx.make_artifact(type='FeatureTable[Zscore]',
                                        view=zscores_filepath,
                                        view_type=PepsirfContingencyTSVFormat)
        
        if negative_data_filepath:
                negative_data = ctx.make_artifact(type='FeatureTable[Normed]',
                                                        view=negative_data_filepath,
                                                        view_type=PepsirfContingencyTSVFormat)
        else:
                negative_data = None

        if highlighted_probes_filepath:
            highlighted_probes = ctx.make_artifact(type='InfoSNPN',
                                                        view=highlighted_probes_filepath,
                                                        view_type=PepsirfInfoSNPNFormat)

        zenrich_vis, = zenrich(data = data,
            zscores = zscores,
            negative_controls = negative_controls,
            highlight_probes = highlighted_probes_filepath,
            source = source,
            pn_filepath = pn_filepath,
            peptide_metadata = peptide_metadata,
            tooltip = tooltip,
            color_by = color_by,
            negative_data = negative_data,
            step_z_thresh = step_z_thresh,
            upper_z_thresh = upper_z_thresh,
            lower_z_thresh = lower_z_thresh,
            exact_z_thresh = exact_z_thresh,
            pepsirf_binary = pepsirf_binary)

        return zenrich_vis
