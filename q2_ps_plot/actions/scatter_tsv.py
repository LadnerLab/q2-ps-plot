from q2_pepsirf.format_types import PepsirfContingencyTSVFormat

def repScatters_tsv(
    ctx,
    source = None,
    pn_filepath= None,
    plot_log= False,
    zscore_filepath = None,
    col_sum_filepath = None,
    facet_charts = False
):
    repScatters = ctx.get_action('ps-plot', 'repScatters')

    # import data into an artifact
    if zscore_filepath:
        zscore = ctx.make_artifact(type='FeatureTable[Zscore]',
                                    view=zscore_filepath,
                                view_type=PepsirfContingencyTSVFormat)
    else:
        zscore = None

    # import data into an artifact
    if col_sum_filepath:
        col_sum = ctx.make_artifact(type='FeatureTable[Normed]',
                                    view=col_sum_filepath,
                                view_type=PepsirfContingencyTSVFormat)
    else:
        col_sum = None

    repScatters_vis, = repScatters(
        source = source,
        pn_filepath = pn_filepath,
        plot_log = plot_log,
        zscore = zscore,
        col_sum = col_sum,
        facet_charts = facet_charts
    )

    return repScatters_vis