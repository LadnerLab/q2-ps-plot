import pandas as pd
import altair as alt
import numpy as np
from collections import defaultdict
import qiime2, itertools, os


# function to collect pairwise pairs of samples
def _make_pairs_list(column):
    series = column.to_series()
    pairs = {k: v.index for k,v in series.groupby(series)}
    result = []
    for _, ids in pairs.items():
        combination = list(itertools.combinations(ids, 2))
        if len(combination) > 0:
            result.append(tuple(combination))
        else:
            result.append(ids)
    return result

def repScatters(
    output_dir: str,
    source: qiime2.CategoricalMetadataColumn,
    plot_log: bool = False,
    zscore: pd.DataFrame = None,
    col_sum: pd.DataFrame = None)->None:

    # check wether zscore matrix or colsum matrix was provided
    if zscore is not None:
        data = zscore.transpose()
    if col_sum is not None:
        data = col_sum.transpose()

    # start the default dict and add the following column names or keys
    hmDic = defaultdict(list)
    columnNms = ['sample', 'bin_x_start', 'bin_x_end', 'bin_y_start', 'bin_y_end', 'count']
    for nm in columnNms:
        hmDic[nm]

    # call function to collect samples pairs
    pairs = _make_pairs_list(source)

    # initialize empty list to collet sample names for dropdown menu
    samples = []

    # start loop to collect heatmap values
    for lst in pairs:
        for tpl in lst:

            #set x and y values from the dataframe
            if type(tpl) == 'tuple':
                samp = '~'.join(tpl)
                samples.append(samp)
                x = list(data[tpl[0]])
                y = list(data[tpl[1]])
            else:
                samp = tpl
                samples.append(samp)
                x = list(data[tpl])
                y = list(data[tpl])

            if plot_log:
                x = np.log10(np.array(x)+1)
                y = np.log10(np.array(y)+1)
                xTitle = "Sample1 log10(score + 1)"
                yTitle = "Sample2 log10(score + 1)"
            else:
                xTitle = "Sample1"
                yTitle = "Sample2"

            #generate heatmap values
            heatmap, xedges, yedges = np.histogram2d(x, y, bins=(70,70))
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

                    #place heatmap values in the dictionary created
                    hmDic['sample'].append(samp)
                    hmDic['bin_x_start'].append(bin_x_start)
                    hmDic['bin_x_end'].append(bin_x_end)
                    hmDic['bin_y_start'].append(bin_y_start)
                    hmDic['bin_y_end'].append(bin_y_end)
                    hmDic['count'].append(count)
    
    #convert the heatmap to a dataframe
    hmDf = pd.DataFrame(hmDic)

    # set the dropdown values
    sample_dropdown = alt.binding_select(options=samples, name='Sample Select')
    sample_select = alt.selection_single(fields=['sample'], bind=sample_dropdown, name="sample", init={'sample': samples[0]})

    # create scatterplot chart with attached dropdown menu
    heatmapChart = alt.Chart(hmDf).mark_rect().encode(
            alt.X('bin_x_start:Q', title=xTitle),
            alt.X2('bin_x_end:Q'),
            alt.Y('bin_y_start:Q', title=yTitle),
            alt.Y2('bin_y_end:Q'),
            alt.Color('count:Q', scale = alt.Scale(scheme='plasma'))
    ).add_selection(
        sample_select
    ).transform_filter(
        sample_select
    )

    # save the chart to index.html to be saved to a .qzv file
    heatmapChart.save(os.path.join(output_dir, "index.html"))


    
