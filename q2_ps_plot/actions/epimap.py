from collections import defaultdict
from q2_pepsirf.format_types import (
    PepsirfInfoSNPNFormat, PeptideFastaFmt, PepsirfContingencyTSVFormat
)
from q2_ps_plot.actions.chart import AltChart

import altair as alt
import glob
import os
import pandas as pd
import scipy
import numpy as np
import random
import math

'''
TODO:
- create format type for input files (later)
    - peptide_seq_filepath: PeptideFastaFmt, zscore_filepath: PepsirfContingencyTSVFormat?

'''


# input all files, enrichment type subset name, and color scheme
def epimap(
    output_dir: str,
    metadata_filepath: str,
    peptide_seq_filepath: str,
    zscore_filepath: str,
    p_thresh: float,
    g1_enrichment_subset: list,
    g2_enrichment_subset: list = None,
    fullname_header: str = "FullName",
    codename_header: str = "CodeName",
    protein_header: str = "Protein",
    category_header: str = "Category",
    alascanpos_header: str = "AlaScanPos",
    include_categories: list = ["Target", "Scaffold", "Adjuvant"],
    horizontal_line_pos: list = [0.01, 0.001],
    z_thresh: float = 0,
    xtick_spacing: int = 20,
    color_by_col: str = "Category",
    color_scheme: str = "dark2",
    enriched_output_dir: str = "epimap-enriched-dir",
    enriched_output_filepath: str = "epimap-peptides.tsv"
    ) -> None:

    # read in metadata
    mF = metadata_filepath
    mdf = pd.read_csv(mF, sep="\t", header=0, index_col=1, keep_default_na=False)

    # read in peptide sequences
    fasta_dict = read_fasta_dict_upper(peptide_seq_filepath)

    # read in zscores
    zF = zscore_filepath
    zdf = pd.read_csv(zF, sep="\t", header=0, index_col=0)
    zdf = zdf.sort_index()

    # get groups
    g1_DF = zdf[[x for x in zdf.columns if any(substr in x for substr in g1_enrichment_subset)]]

    if g2_enrichment_subset == None:
        g2_DF = zdf[[x for x in zdf.columns if all(substr not in x for substr in g1_enrichment_subset)]]
    else:
        g2_DF = zdf[[x for x in zdf.columns if any(substr in x for substr in g2_enrichment_subset)]]


    # Collect info for a volcano plot
    thresh = 10
    minHits = 2

    chart_dict = defaultdict(list)

    for pn, Zs in g1_DF.iterrows():
        #Want to exclude certain categories of peptides, to focus in on the normal sliding window design
        if mdf[category_header][pn] in include_categories and not mdf[alascanpos_header][pn]:
            hits = sum([x>thresh for x in Zs])
            if hits >= minHits:
                res = scipy.stats.ttest_ind(Zs, g2_DF.loc[pn], equal_var=False, alternative="greater")

                chart_dict['CodeName'].append(pn)
                chart_dict['FullName'].append(mdf[fullname_header][pn])
                chart_dict['Protein'].append(mdf[protein_header][pn])
                chart_dict['Category'].append(mdf[category_header][pn])
                chart_dict['Seq'].append(fasta_dict[pn])
                chart_dict['ZScoreDiff'].append(np.average(Zs) - np.average(g2_DF.loc[pn]))
                chart_dict['PVal'].append(res.pvalue)

    chartDf = pd.DataFrame.from_dict(chart_dict)

    chart = AltChart(
        dataframe=chartDf,
        mark_size=60,
        x_val="ZScoreDiff:Q",
        x_title="Zscore Difference",
        x_axis_ticks=[z_thresh] + list(range(
                        round("down", (chartDf["ZScoreDiff"].min())), 
                        round("up", (chartDf["ZScoreDiff"].max())), 
                        xtick_spacing
                        )),
        y_val="PVal:Q",
        y_title="T-test p-value",
        y_axis_ticks=[p_thresh] + horizontal_line_pos,
        color_by=f"{color_by_col}:N",
        color_scheme=color_scheme,
        tooltip=[
            "CodeName:N",
            "FullName:N",
            "Protein:N",
            "Category:N",
            "Seq:N"
            ]
        ).make_epimap()
    

    y_line = alt.Chart(pd.DataFrame({'y': [p_thresh]})).mark_rule(strokeDash=[2,1], strokeWidth=1).encode(y='y')
    x_line = alt.Chart(pd.DataFrame({'x': [z_thresh]})).mark_rule(strokeDash=[2,1], strokeWidth=1).encode(x='x')

    final_plot = chart + x_line + y_line

    final_plot = final_plot.configure_view(
                                continuousHeight=500,
                                continuousWidth=500
                                )

    final_plot.save(os.path.join(output_dir, "index.html"))

    # save peptide data to outfile (tsv)
    threshDf = chartDf.loc[(chartDf['ZScoreDiff'] >= z_thresh) & (chartDf['PVal'] <= p_thresh)]
    threshDf[["CodeName", "ZScoreDiff", "PVal"]].to_csv(enriched_output_filepath, sep="\t", index=False)

    # setup enriched dir file to use with heatmap
    if not os.path.exists(enriched_output_dir):
        os.mkdir(enriched_output_dir)

    open( f"{enriched_output_dir}/failedEnrichment.txt", "w" ).close()
    # save the peptide names to directory file
    threshDf["CodeName"].to_csv(f"{enriched_output_dir}/subset_enriched.txt", index=False, header=False)


# source: https://github.com/jtladner/Modules
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

# source: https://github.com/jtladner/Modules
def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

def round(type: str, x: float) -> int:
    if type=="up":
        return math.ceil(x / 10.0) * 10
    if type=="down":
        return math.floor(x / 10.0) * 10





