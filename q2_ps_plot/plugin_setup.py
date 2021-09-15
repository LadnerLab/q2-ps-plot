#!/usr/bin/env python
import importlib
import q2_ps_plot

from qiime2.plugin import (Plugin,
                        SemanticType,
                        model,
                        Int,
                        Range,
                        MetadataColumn,
                        Categorical,
                        Str,
                        List,
                        Visualization,
                        Metadata)

from q2_ps_plot.format_types import (
    Normed, Zscore, PepsirfContingencyTSVDirFmt, PepsirfContingencyTSVFormat)
import q2_ps_plot.actions as actions

from q2_types.feature_table import FeatureTable, BIOMV210DirFmt


# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin("ps-plot", version=q2_ps_plot.__version__,
                website="https://github.com/LadnerLab/q2-ps-plot")

plugin.register_formats(PepsirfContingencyTSVFormat,
                        PepsirfContingencyTSVDirFmt)

plugin.register_semantic_types(Normed, Zscore)
plugin.register_semantic_type_to_format(FeatureTable[Normed | Zscore],
                                        BIOMV210DirFmt)

shared_parameters = {
    "step_z_thresh": Int % Range(1, None),
    "upper_z_thresh": Int % Range(2, None),
    "lower_z_thresh": Int % Range(1, None),
    "exact_z_thresh": List[Str],
    "source": MetadataColumn[Categorical],
    "peptide_metadata": Metadata,
    "pepsirf_binary": Str,
    "negative_controls": List[Str],
    "tooltip": List[Str]
}
shared_descriptions = {
    "step_z_thresh": "Integar to increment z-score thresholds.",
    "upper_z_thresh": "Upper limit of z-score thresholds (non-inclusive).",
    "lower_z_thresh": "Lower limit of z-score thresholds (inclusive).",
    "exact_z_thresh": "List of exact z score thresholds either individual or combined. "
                    "List MUST BE in descending order. (Example argument: '--p-exact-z-thresh 25 10 3' "
                    "or '--p-exact-z-thresh 6,25 4,10 1,3')",
    "source": "Metadata file containing all sample names and their source groups. "
            "Used to create pairs tsv to run pepsirf enrich module.",
    "peptide_metadata": "Filename of file that contains peptide metadata related to "
                        "data to be plotted.",
    "pepsirf_binary": "The binary to call pepsirf on your system. "
                    "Used to call pepsirf enrich module.",
    "negative_controls": "Sample names of the negative controls to be used "
                        "(Example argument: --p-negative-controls sample1 sample2 sample3).",
    "tooltip": "List of title names found in the peptide metadata file "
                "to be added to the hover tooltip (Parameter is case sensitive. "
                "Example argument: --p-tooltip Species SpeciesID). "
                "'Peptide' and 'Zscores' will always be added to the list of titles, "
                "if peptide metadata is not provided just 'Peptide' and 'Zscores' will be shown."
}

plugin.pipelines.register_function(
    function=actions.zenrich_tsv,
    inputs={},
    outputs=[
        ("zenrich_vis", Visualization)
    ],
    parameters={
        'data_filepath': Str,
        'zscores_filepath': Str,
        'negative_data_filepath': Str,
        **shared_parameters
    },
    input_descriptions=None,
    output_descriptions=None,
    parameter_descriptions={
        'data_filepath': "Filepath of .tsv file containing normalized read counts of samples and peptides. "
                    "First column header must be 'Sequence Name' as produced by pepsirf.",
        'zscores_filepath': "Filepath of .tsv file containing z scores of the normalized read counts. "
                    "Fist column header must be 'Sequence Name' as produced by pepsirf.",
        'negative_data_filepath':"Filepath of .tsv file containing normalized read counts of controls and peptides. "
                            "First column header must be 'Sequence Name' as produced by pepsirf.",
        **shared_descriptions
    },
    name='zenrich TSV Pipeline',
    description="Pipeline that converts .tsv files to .qza files and then runs zenrich."
)

plugin.visualizers.register_function(
    function=actions.zenrich,
    inputs={
        'data': FeatureTable[Normed],
        'zscores': FeatureTable[Zscore],
        'negative_data': FeatureTable[Normed]
    },
    parameters=shared_parameters,
    input_descriptions={
        'data': "FeatureTable containing normalized read counts of samples and peptides. "
                "First column header must be 'Sequence Name' as produced by pepsirf.",
        'zscores': "FeatureTable containing z scores of the normalized read counts. "
                "Fist column header must be 'Sequence Name' as produced by pepsirf.",
        'negative_data': "FeatureTable containing normalized read counts of controls and peptides. "
                        "First column header must be 'Sequence Name' as produced by pepsirf."
    },
    parameter_descriptions=shared_descriptions,
    name='Z Enrichment Variance Visualizer',
    description="Creates a scatterplot of enriched peptides, points are colored "
                "according to the z score thresholds provided. Scatterplot is "
                "layered over a heatmap containing all of the data."
)



importlib.import_module("q2_ps_plot.transformers")
