#!/usr/bin/env python
import pandas as pd
from q2_ps_plot.format_types import PepsirfContingencyTSVFormat
import qiime2

def zenrich_tsv(ctx,
            data_filepath,
            zscores_filepath,
            negative_controls,
            source = None,
            pn_filepath = None,
            peptide_metadata = None,
            tooltip = ['Species', 'SpeciesID'],
            negative_data_filepath = None,
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

        zenrich_vis, = zenrich(data = data,
            zscores = zscores,
            negative_controls = negative_controls,
            source = source,
            pn_filepath = pn_filepath,
            peptide_metadata = peptide_metadata,
            tooltip = tooltip,
            negative_data = negative_data,
            step_z_thresh = step_z_thresh,
            upper_z_thresh = upper_z_thresh,
            lower_z_thresh = lower_z_thresh,
            exact_z_thresh = exact_z_thresh,
            pepsirf_binary = pepsirf_binary)

        return zenrich_vis
