import pandas as pd
import altair as alt
import os
from altair_saver import save

from q2_pepsirf.format_types import PepsirfInfoSumOfProbesFmt

def readCountsBoxplot(
    output_dir: str,
    read_counts: pd.DataFrame) -> None:
    
    #collects min and max values of sum of probes
    rc_min = min(read_counts['Sum of probe scores'])
    rc_max = max(read_counts['Sum of probe scores'])


    #creates boxplot
    chart = alt.Chart(read_counts).transform_fold(
        ['Sum of probe scores'],
        as_=['key', 'value']
    ).mark_boxplot(size = 50).encode(
        x=alt.X('key:N', axis = alt.Axis(title="Sum Of Probes")),
        y=alt.Y('value:Q', axis = alt.Axis(title="Value (1e6)"), scale=alt.Scale(domain=[rc_min, rc_max]))
    ).properties(width=300)

    # saves boxplot as a png file and to index.html for creation of qzv file
    save(chart, "readCountBoxplot.png", scale_factor=10)
    chart.save(os.path.join(output_dir, "index.html"))

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

    # create chart for boxplot
    chart = alt.Chart(sumDf).transform_fold(
        ['sum of enriched'],
        as_=['key', 'value']
    ).mark_boxplot(size = 50).encode(
        x=alt.X('key:N', axis = alt.Axis(title="Sum Of Enriched")),
        y=alt.Y('value:Q', axis = alt.Axis(title="Value"))
    ).properties(width=300)

    # saves boxplot as a png file and to index.html for creation of qzv file
    save(chart, "enrichedRCBoxplot.png", scale_factor=10)
    chart.save(os.path.join(output_dir, "index.html"))