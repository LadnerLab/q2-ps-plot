import pandas as pd
import altair as alt
import os
from altair_saver import save

# function that creates and saves the boxplot as a png and qzv
def _make_box_plot(dataframe, key, output_dir, png_out, rc_min, rc_max):
    #creates boxplot
    chart = alt.Chart(dataframe).transform_fold(
        [key],
        as_=['key', 'value']
    ).mark_boxplot(size = 50).encode(
        x=alt.X('key:N', axis = alt.Axis(title="Sum Of Probes")),
        y=alt.Y('value:Q', axis = alt.Axis(title="Value"), scale=alt.Scale(domain=[rc_min, rc_max]))
    ).properties(width=300)

    # saves boxplot as a png file and to index.html for creation of qzv file
    save(chart, png_out, scale_factor=10)
    chart.save(os.path.join(output_dir, "index.html"))

def readCountsBoxplot(
    output_dir: str,
    read_counts: pd.DataFrame) -> None:
    
    #collects min and max values of sum of probes for scaling of graph
    rc_min = min(read_counts['Sum of probe scores'])
    rc_max = max(read_counts['Sum of probe scores'])

    # create and save boxplot
    _make_box_plot(read_counts, 'Sum of probe scores', output_dir, "readCountBoxplot.png", rc_min, rc_max)

def enrichmentRCBoxplot(
    output_dir: str,
    enriched_dir: pd.DataFrame) -> None:

    # collect file names
    files = list(enriched_dir.columns)

    # create dictionary for collection of sum of enriched
    sumDict = {"sum of enriched": []}

    # loop through the file names and collect sum of True values for enriched
    for f in files:
        sumDict["sum of enriched"].append(len(enriched_dir[enriched_dir[f] == True]))

    # convert dictionary to dataframe
    sumDf = pd.DataFrame(sumDict)

    #find max/min for scaleing of graph
    rc_min = min(sumDf['sum of enriched'])
    rc_max = max(sumDf['sum of enriched'])

    # create and save boxplot
    _make_box_plot(sumDf, 'sum of enriched', output_dir, "enrichedCountBoxplot.png", rc_min, rc_max)