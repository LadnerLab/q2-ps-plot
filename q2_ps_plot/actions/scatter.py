import pandas as pd
import altair as alt
import numpy as np
from collections import defaultdict
import qiime2, itertools, os

def _make_pairs_list(column):
    series = column.to_series()
    pairs = {k: v.index for k,v in series.groupby(series)}
    result = []
    for _, ids in pairs.items():
        result.append(list(itertools.combinations(ids, 2)))
    return result

def repZscatters(
    output_dir: str,
    zscore: pd.DataFrame,
    source: qiime2.CategoricalMetadataColumn = None)->None:

    zscore = zscore.transpose()

    hmDic = defaultdict(list)
    columnNms = ['sample', 'bin_x_start', 'bin_x_end', 'bin_y_start', 'bin_y_end', 'count']
    for nm in columnNms:
        hmDic[nm]

    pairs = _make_pairs_list(source)
    samples = []
    for lst in pairs:
        for tpl in lst:
            samp = '~'.join(tpl)
            samples.append(samp)
            x = list(zscore[tpl[0]])
            y = list(zscore[tpl[1]])

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

                    hmDic['sample'].append(samp)
                    hmDic['bin_x_start'].append(bin_x_start)
                    hmDic['bin_x_end'].append(bin_x_end)
                    hmDic['bin_y_start'].append(bin_y_start)
                    hmDic['bin_y_end'].append(bin_y_end)
                    hmDic['count'].append(count)
    
    hmDf = pd.DataFrame(hmDic)

    sample_dropdown = alt.binding_select(options=samples, name='Sample Select')
    sample_select = alt.selection_single(fields=['sample'], bind=sample_dropdown, name="sample", init={'sample': samples[0]})

    heatmapChart = alt.Chart(hmDf).mark_rect().encode(
            alt.X('bin_x_start:Q', title="Sample1"),
            alt.X2('bin_x_end:Q'),
            alt.Y('bin_y_start:Q', title="Sample2"),
            alt.Y2('bin_y_end:Q'),
            alt.Color('count:Q', scale = alt.Scale(scheme='plasma'))
    ).add_selection(
        sample_select
    ).transform_filter(
        sample_select
    )

    heatmapChart.save(os.path.join(output_dir, "index.html"))


    
