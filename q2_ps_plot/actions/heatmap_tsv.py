from q2_pepsirf.format_types import EnrichedPeptideDirFmt, ProteinAlignmentFmt

def proteinHeatmap_dir(
        ctx,
        enriched_dir_filepath,
        protein_alignment_filepath,
        enriched_suffix="_enriched.txt",
        align_header="AlignPos",
        align_delim="~",
        include_species=False,
        species_header="Species",
        color_scheme="viridis",
        output_size=500
        ):
    
    proteinHeatmap = ctx.get_action("ps-plot", "proteinHeatmap")
    
    enriched_dir = ctx.make_artifact(
        type="PairwiseEnrichment",
        view=enriched_dir_filepath,
        view_type=EnrichedPeptideDirFmt
    )

    protein_alignment = ctx.make_artifact(
        type="ProteinAlignment",
        view=protein_alignment_filepath,
        view_type=ProteinAlignmentFmt
    )

    proteinHeatmap_vis, = proteinHeatmap(
        enriched_dir=enriched_dir,
        protein_alignment=protein_alignment,
        enriched_suffix=enriched_suffix,
        align_header=align_header,
        align_delim=align_delim,
        include_species=include_species,
        species_header=species_header,
        color_scheme=color_scheme,
        output_size=output_size
    )

    return proteinHeatmap_vis

