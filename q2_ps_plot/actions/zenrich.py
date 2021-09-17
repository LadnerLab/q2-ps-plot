#!/usr/bin/env python
import pandas as pd
import numpy as np
import altair as alt
import glob, os, subprocess
import tempfile
import csv
from collections import defaultdict

from q2_ps_plot.format_types import PepsirfContingencyTSVFormat, Zscore
import qiime2
from q2_types.feature_table import BIOMV210Format

def _make_pairs_file(column, outpath):
    series = column.to_series()
    pairs = {k: v.index for k,v in series.groupby(series)}
    with open(outpath, 'w') as fh:
        for _, ids in pairs.items():
            fh.write('\t'.join(ids) + '\n')

def zenrich(output_dir: str,
            data: pd.DataFrame,
            zscores: PepsirfContingencyTSVFormat,
            negative_controls: list,
            source: qiime2.CategoricalMetadataColumn = None,
            pn_filepath: str = None,
            peptide_metadata: qiime2.Metadata = None,
            tooltip: list=['Species', 'SpeciesID'],
            negative_data: pd.DataFrame = None,
            step_z_thresh: int=5,
            upper_z_thresh: int=30,
            lower_z_thresh: int=5,
            exact_z_thresh: list=None,
            pepsirf_binary: str="pepsirf") -> None:

    old = os.getcwd() # TODO: bug in framework, remove this when fixed

    #create temporary directory to work in
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)

        #create pairs file
        if source:
            pairsFile = os.path.join(tempdir, 'pairs.tsv')
            _make_pairs_file(source, pairsFile)
        elif pn_filepath:
            pairsFile = os.path.abspath(pn_filepath)

        #flip data frame
        data = data.transpose()

        #put zscores data into dataframe and flip
        zBiom = zscores.view(BIOMV210Format)
        zData = zBiom.view(pd.DataFrame)
        zData = zData.transpose()

        #set negative matrix
        if not negative_data:
            negative_data = data
        else:
            negative_data = negative_data.transpose()

        #set max rows for altair to none
        alt.data_transformers.enable('default', max_rows=None)
        
        #values created by user input
        threshFile = "tempThreshFile.tsv"
        
        #set values for pepsirf
        outSuffix = "_tempEnriched.txt"
        outputDir = []
        
        #convert colsum data to pandas data frame
        peptideNames = negative_data.index
        
        #get list of thresholds
        if not exact_z_thresh:
            threshRange = list(reversed(range(lower_z_thresh, upper_z_thresh, step_z_thresh)))
        else:
            threshRange = exact_z_thresh

        #iterate through thresholds and generate threshold file to be used to get enriched peptide
        for score in threshRange:
            outputDir.append(str(score))
            #write the threshold file for specified zscore
            with open(threshFile, 'w', newline='') as out_file:
                    tsv_writer = csv.writer(out_file, delimiter='\t')
                    tsv_writer.writerow([str(zscores), score])
            
            #run p enrich module
            cmd = "%s enrich -t %s -s %s -x %s -o %s >> enrich.out" % (pepsirf_binary, threshFile, pairsFile, outSuffix, str(score))
            subprocess.run(cmd, shell=True)

        #create empty dicionaries to collect all sample names and sample/peptide combos
        dropNames = {}
        pepSamp = {}

        #create empty dictionary of lists to collect enriched peptides
        enrDict = defaultdict(list)
        columnNms = ['Peptide','sample','sample_value', 'negative_control', 'z_score_threshold', 'Zscores']
        for nm in columnNms:
            enrDict[nm]

        #iterate through created directories
        for oD in outputDir:
            enrFiles = glob.glob("%s/*%s" % (oD, outSuffix))
            
            #iterate through enriched files to gather the peptides
            for file in enrFiles:
                
                #get individual sample names
                sNames = os.path.basename(file).split(outSuffix)[0].split("~")
                
                #get joined sample name
                sample = "~".join(sNames)

                #generate a list of sample names
                if sample not in dropNames:
                    dropNames[sample] = ""
                
                #collect information for the pandas data frame
                with open(file, 'r') as fin:
                    for row in fin:
                        pep = row.rstrip("\n")
                        if (pep,sample) not in pepSamp:
                            x = np.mean([float(negative_data[sn][pep]) for sn in negative_controls])
                            y = np.mean([float(data[sn][pep]) for sn in sNames])
                            z = [zData[sn][pep] for sn in sNames]
                            zToStr = ', '.join([str(elem) for elem in z])

                            #add all enriched peptide info to enriched dictionary
                            enrDict['Peptide'].append(pep)
                            enrDict['sample'].append(sample)
                            enrDict['sample_value'].append(np.log10(y+1))
                            enrDict['negative_control'].append(np.log10(x+1))
                            enrDict['z_score_threshold'].append(oD)
                            enrDict['Zscores'].append(zToStr)
                            pepSamp[(pep,sample)] = ""

        #convert eriched dictionary to enriched data frame
        enrichedDf = pd.DataFrame(enrDict)

        #collect peptide metadata information and add it to the enriched data frame
        if peptide_metadata:
            metaDf = peptide_metadata.to_dataframe()
            enrichedDf = metaDf.merge(enrichedDf, how="right", right_on="Peptide", left_index=True)

            #add pepsirf to tooltip list and add data type to each item
            tooltip.insert(0, 'Peptide')
            tooltip.insert(1, 'Zscores')
            i = 0
            for title in tooltip:
                tooltip[i] += ':N'
                i += 1
        else:
            tooltip = ['Peptide:N', 'Zscores:N']

        #create empty directory to collect all heatmap information
        hmDic = defaultdict(list)
        columnNms = ['sample', 'bin_x_start', 'bin_x_end', 'bin_y_start', 'bin_y_end', 'count']
        for nm in columnNms:
            hmDic[nm]

        #collect x axis data
        x = np.array([np.mean([float(negative_data[sn][pn]) for sn in negative_controls]) for pn in peptideNames])
        xLog = np.log10(x+1)

        #iterate through samples to fill Data Frame
        for sample in dropNames.keys():
            sNames = sample.split('~')
            y = np.array([np.mean([float(data[sn][pn]) for sn in sNames]) for pn in peptideNames])
            yLog = np.log10(y+1)
            heatmap, xedges, yedges = np.histogram2d(xLog, yLog, bins=(70,70))

            #fill binned data frame for heatmap generation by row
            for x in range(0, heatmap.shape[0]):
                for y in range(0, heatmap.shape[1]):
                    count = heatmap[x,y]
                    #do not add data with count of 0
                    if count == 0.0:
                        continue
                    bin_x_start = xedges[x]
                    bin_x_end = xedges[x+1]
                    bin_y_start = yedges[y]
                    bin_y_end = yedges[y+1]

                    #add all x and y bins to heatmap dictionary
                    hmDic['sample'].append(sample)
                    hmDic['bin_x_start'].append(bin_x_start)
                    hmDic['bin_x_end'].append(bin_x_end)
                    hmDic['bin_y_start'].append(bin_y_start)
                    hmDic['bin_y_end'].append(bin_y_end)
                    hmDic['count'].append(count)

        #convert heatmap dictionary to data frame
        heatmapDf = pd.DataFrame(hmDic)

        #find axis ratio for chart width an height
        xy_max = heatmapDf[['bin_x_end', 'bin_y_end']].max()
        ratio = (xy_max[0] / xy_max[1]) - 1 
        chartHeight = 500
        chartWidth =  chartHeight + (20 * ratio)

        #create dropdown specs
        sample_dropdown = alt.binding_select(options=list(dropNames.keys()), name='Sample Select')
        sample_select = alt.selection_single(fields=['sample'], bind=sample_dropdown, name="sample", init={'sample': list(dropNames.keys())[0]})
        
        #create scatterplot of enriched peptides
        scatter = alt.Chart(enrichedDf).mark_circle(size=50).encode(
            x = alt.X('negative_control:Q', title = "Negative Control log10(value+1.0)"),
            y = alt.Y('sample_value:Q', title = "Sample log10(value+1.0)"),
            color = alt.Color('z_score_threshold:N',
                scale=alt.Scale(range=['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']),
                sort=threshRange,
                legend=alt.Legend(title='Z Score Thresholds')),
            tooltip = tooltip
        ).add_selection(
            sample_select
        ).transform_filter(
            sample_select
        )

        #create heatmap of all peptides
        heatmapChart = alt.Chart(heatmapDf, width=chartWidth, height=chartHeight).mark_rect().encode(
            alt.X('bin_x_start:Q', scale=alt.Scale(domain=(0,xy_max[0]))),
            alt.X2('bin_x_end:Q'),
            alt.Y('bin_y_start:Q', scale=alt.Scale(domain=(0,xy_max[1]))),
            alt.Y2('bin_y_end:Q'),
            alt.Color('count:Q', scale = alt.Scale(scheme='greys'), legend=alt.Legend(title='Point Frequency'))
        ).transform_filter(
            sample_select
        )

        #layer scatterplot and heatmap
        finalChart = alt.layer(heatmapChart, scatter).properties(title="Z Score Threshold Variance")

        #save the plot
        finalChart.save(os.path.join(output_dir, "index.html"))

    os.chdir(old)  # TODO: found a bug in the framework, remove when fixed