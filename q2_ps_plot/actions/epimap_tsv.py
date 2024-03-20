#!/usr/bin/env python
from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSNPNFormat
)

import qiime2


def epimap_dir(
    ctx,
    metadata_filepath,
    peptide_seq_filepath,
    zscore_filepath,
    p_thresh,
    g1_enrichment_subset,
    g2_enrichment_subset=None,
    fullname_header="FullName",
    codename_header="CodeName",
    protein_header="Protein",
    category_header="Category",
    alascanpos_header="AlaScanPos",
    include_categories=["Target", "Scaffold", "Adjuvant"],
    horizontal_line_pos=[0.01, 0.001],
    color_by_col="Category",
    color_scheme="dark2"
    ):

    # collect the epimpap module action
    epimap = ctx.get_action("ps-plot", "epimap")


    # run the epimap module with all the inputs and parameters given
    epimap_vis, = epimap(
        metadata_filepath=metadata_filepath,
        peptide_seq_filepath=peptide_seq_filepath,
        zscore_filepath=zscore_filepath,
        p_thresh=p_thresh,
        g1_enrichment_subset=g1_enrichment_subset,
        g2_enrichment_subset=g2_enrichment_subset,
        fullname_header=fullname_header,
        codename_header=codename_header,
        protein_header=protein_header,
        category_header=category_header,
        alascanpos_header=alascanpos_header,
        include_categories=include_categories,
        horizontal_line_pos=horizontal_line_pos,
        color_by_col=color_by_col,
        color_scheme=color_scheme
    )

    #return the visualization
    return epimap_vis
