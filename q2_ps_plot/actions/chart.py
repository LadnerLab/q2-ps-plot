import altair as alt
import pandas as pd

# TODO: make this a class

# Name: make_alt_chart
# Process: creates a chart given information on what kind of chart to create
#			and what information to fill it
# Method inputs/parameters: 
# Method outputs/returned: altair chart
# Dependencies: pandas, altair
def make_alt_chart( 
	chart_type: str,
	dataframe: pd.DataFrame,
	x_val: str, 
	x_title: str, 
	y_val: str, 
	y_title: str,
	tooltip: list,
	x_max: int = 6,
	x_axis_vals: list = None,
	y_axis_vals: list = None,
	y_scale_setup: alt.Scale = None,
	color_by: str = None, 
	color_scheme: str = None,
	color_range: list = None
	):
	# default values
	x_bin_data = None
	x_axis_data = None
	y_axis_data = None

	chart = alt.Chart(dataframe)

	# chart types
	if chart_type == "heatmap":
		chart = chart.mark_rect()
		x_bin_data = alt.Bin(
	            		maxbins=x_max,
	            		minstep=1
            		)

	elif chart_type == "scatterplot":
		chart = chart.mark_circle(size=50)

		if x_axis_vals != None:
			x_axis_data = alt.Axis(
						values = x_axis_vals
						)
		if y_axis_vals != None:
			y_axis_data = alt.Axis(
						values = y_axis_vals
						)



	# encode information
	chart = chart.encode(
        x=alt.X(
            x_val,
            title=x_title,
            bin=x_bin_data,
            scale=alt.Scale(
            	zero=True
            	),
            axis=x_axis_data
        ),

        y=alt.Y(
            y_val,
            title=y_title,
            axis=y_axis_data,
            scale=y_scale_setup
        ),
        color=alt.Color(
            color_by,
            scale=alt.Scale(
                scheme=color_scheme
            )
        ),
        tooltip=tooltip
    )

	return chart