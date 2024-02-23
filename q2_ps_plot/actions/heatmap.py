from collections import defaultdict
from q2_pepsirf.format_types import EnrichedPeptideDirFmt, ProteinAlignmentDirFormat

import altair as alt
import glob
import os
import pandas as pd


# Name: proteinHeatmap
# Process: creates an interactive heatmap to show the alignment of proteins
# Method inputs/parameters: output_dir, enriched_dir, protein_alignment,
# enriched_suffix, start_header, stop_header
# Dependencies: glob, os, pandas, altair, and defaultdict
def proteinHeatmap(
        output_dir: str,
        enriched_dir: EnrichedPeptideDirFmt,
        protein_alignment: ProteinAlignmentDirFormat,
        enriched_suffix: str = "_enriched.txt",
        align_header: str = "AlignPos",
        align_delim: str = "~",
        include_species: bool = False,
        species_header: str = "Species",
        color_scheme: str = "viridis") -> None:

    print(os.listdir(str(protein_alignment)))   

    # collect enriched peptide files
    enrFiles = glob.glob("%s/*%s" % (str(enriched_dir), enriched_suffix))

    # create a dictionary with list values to collect peptide information
    enrDict = defaultdict(list)

    # go through each enriched peptide file and collect peptide information
    for file in enrFiles:
        sName = os.path.basename(file).split(enriched_suffix)[0]
        with open(file, "r") as fin:
            for row in fin:
                pep = row.rstrip("\n")
                enrDict[sName].append(pep)

    # put the enriched peptide directory into a pandas DataFrame
    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in enrDict.items()]))
    df = df.melt(var_name="sample", value_name="peptide")

    # create a dictionary of a dictionary to collect alignment/protein/species information
    alignDict = defaultdict()
    proteinDict = defaultdict(dict)

    if include_species:
        speciesDict = defaultdict(dict)

    # open the alignment files and collect the protein, and alignment information, also species
    with open(str(protein_alignment.path /"manifest.tsv"), "r") as fin:
        ct = 0
        for row in fin:
            if ct != 0:
                line = row.rstrip("\n")
                protein, file = line.split("\t")

                # put the alignment file into a pandas DataFrame to easily
                # collect alignment information
                proteinDf = pd.read_csv(
                    str(protein_alignment.path /file),
                    sep="\t", index_col=0
                )

                align = dict(proteinDf[align_header])

                # collect the alignment and the protein
                for peptide, alignment in align.items():
                    alignDict[peptide] = alignment.split(align_delim)
                    proteinDict[peptide]["protein"] = protein

                if include_species:
                    # test if file contains species column
                    if not species_header in proteinDf.columns:
                        raise ValueError(f"\'{species_header}\' column not found in alignment file: \'{file}\'")

                    species = dict(proteinDf[species_header])

                    # collect species
                    for peptide, species in species.items():
                        speciesDict[peptide] = species;

            ct += 1

    # put the alignment dictionary into a pandas dataframe and merge with the
    # enriched dataframe
    alignDf = pd.DataFrame(
        dict([(k,pd.Series(v)) for k,v in alignDict.items()])
    )

    alignDf = alignDf.melt(var_name="peptide", value_name="x")
    newDf = df.merge(
        alignDf, how="right",
        right_on="peptide", left_on="peptide"
    )

    # drop rows with no x value
    newDf.dropna(subset=["x"], inplace=True)

    # put the protein dictionary into a pandas dataframe and merge with the
    # new dataframe
    proteinDf = pd.DataFrame(proteinDict)
    proteinDf = proteinDf.transpose().reset_index()
    proteinDf = proteinDf.rename(columns={"index": "peptide"})
    newDf = newDf.merge(
        proteinDf, how="right",
        right_on="peptide", left_on="peptide"
    )

    if include_species:
        # put the species dictionary into pandas dataframe and merge with newDf
        speciesDf = pd.DataFrame(speciesDict, index=["species"])
        speciesDf = speciesDf.transpose().reset_index()
        speciesDf = speciesDf.rename(columns={"index": "peptide"})
        newDf = newDf.merge(
            speciesDf, how="right", 
            right_on="peptide", left_on="peptide"
        )

    # remove rows with no sample
    newDf.dropna(subset=["sample"], inplace=True)

    count_grouping=["x", "protein", "sample"]
    label_grouping=["x", "protein", "sampleDummy"]

    # remove species if not specified
    if include_species:
        count_grouping.append("species")
        label_grouping.append("species")

    #Create a new sample column to handle NaN values in "sample"
    newDf["sampleDummy"] = newDf["sample"].fillna("dummy")
    # create a count of the x values
    newDf["count"] = newDf.groupby(
        count_grouping
    )["x"].transform("count")
    # Create new label column with all of the peptides for each protein,
    # sample and position combined
    newDf["label"] = newDf.groupby(
        label_grouping
    )["peptide"].transform(lambda x: ",".join(x))
    # Get rid of the 'peptide' column and then remove duplicate rows
    newDf = newDf.drop(columns=["peptide", "sampleDummy"]).drop_duplicates()
    
    # convert type of x column into int and collect the max x value
    newDf["x"] = newDf["x"].astype(int)
    x_max = max(newDf["x"])

    # collect all of the protein values
    proteinList = list(newDf.protein.unique())

    # create the dropdown menus for the chart
    protein_dropdown = alt.binding_select(
        options=proteinList, name="Protein Select"
    )
    protein_select = alt.selection_point(
        fields=["protein"], bind=protein_dropdown,
        name="protein", value=[{"protein": proteinList[0]}]
    )

    if include_species:
        # collect all sample values
        sampleList = list(newDf["sample"].unique())

        # create sample dropdown menu
        sample_dropdown = alt.binding_select(
            options=sampleList, name="Sample Select"
        )
        sample_select = alt.selection_point(
            fields=["sample"], bind=sample_dropdown,
            name="sample", value=[{"sample": sampleList[0]}]
        )

    #set max rows for altair to none
    alt.data_transformers.enable("default", max_rows=None)

    if include_species:
        # generate the heatmap with species as y
        chart = alt.Chart(newDf).mark_rect().encode(
            alt.X(
                "x:Q",
                title="Alignment",
                bin=alt.Bin(maxbins=x_max, minstep=1),
                scale=alt.Scale(zero=True)
            ),

            alt.Y(
                "species:N",
                title="Species"
            ),
            alt.Color(
                "count:Q",
                scale=alt.Scale(
                    scheme=color_scheme
                )
            ),
            tooltip=["label:N"]
        ).add_params(
            sample_select
        ).transform_filter(
            sample_select
        ).add_params(
            protein_select
        ).transform_filter(
            protein_select
        )
    else:
        # generate heatmap with sample as y
        chart = alt.Chart(newDf).mark_rect().encode(
            alt.X(
                "x:Q",
                title="Alignment",
                bin=alt.Bin(maxbins=x_max, minstep=1),
                scale=alt.Scale(zero=True)
            ),
            alt.Y(
                "sample:N",
                title="Samples"
            ),
            alt.Color(
                "count:Q",
                scale=alt.Scale(
                    scheme=color_scheme
                )
            ),
            tooltip=["label:N"]
        ).add_params(
            protein_select
        ).transform_filter(
            protein_select
        )

    # save the chart into the qzv file index
    chart.save(os.path.join(output_dir, "index.html"), scale_factor=10.0)

