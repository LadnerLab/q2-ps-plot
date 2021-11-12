import pandas as pd
import altair as alt
import os
from altair_saver import save

from q2_pepsirf.format_types import PepsirfInfoSumOfProbesFmt

def readCountsBoxplot(
    output_dir: str,
    read_counts: pd.DataFrame) -> None:
    
    # read_counts['all samples'] = read_counts['Sum of probe scores'] / 1000000

    rc_min = min(read_counts['Sum of probe scores'])
    rc_max = max(read_counts['Sum of probe scores'])

    chart = alt.Chart(read_counts).transform_fold(
        ['Sum of probe scores'],
        as_=['key', 'value']
    ).mark_boxplot(size = 50).encode(
        x=alt.X('key:N', axis = alt.Axis(title="Sum Of Probes")),
        y=alt.Y('value:Q', axis = alt.Axis(title="Value (1e6)"), scale=alt.Scale(domain=[rc_min, rc_max]))
    ).properties(width=300)

    save(chart, "readCountBoxplot.png")
    chart.save(os.path.join(output_dir, "index.html"))