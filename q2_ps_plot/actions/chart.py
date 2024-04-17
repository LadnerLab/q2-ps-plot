import altair as alt
import pandas as pd

# store the various chart creations here for reuse
class AltChart:

	def __init__(
		self,
		dataframe: pd.DataFrame,
		x_val: str, 
		x_title: str, 
		y_val: str, 
		y_title: str,
		tooltip: list,
		x2_val: str = None,
		y2_val: str = None,
		x_axis_ticks: list = None,
		y_axis_ticks: list = None,
		x_max: int = 6,
		mark_size: int = None,
		legend_title = None,
		color_by: str = None, 
		color_scheme: str = None,
		color_range: list = None
	):
		self.df = dataframe
		self.x_val = x_val
		self.x_title = x_title
		self.y_val = y_val
		self.y_title = y_title
		self.tooltip = tooltip
		self.x2_val = x2_val
		self.y2_val = y2_val
		self.x_axis_ticks = x_axis_ticks
		self.y_axis_ticks = y_axis_ticks
		self.x_max = x_max
		self.mark_size = mark_size
		self.legend_title = legend_title
		self.color_by = color_by
		self.color_scheme = color_scheme
		self.color_range = color_range

	# create chart for heatmap plugin
	def make_heatmap(self):
		chart = alt.Chart(self.df).mark_rect().encode(
			alt.X(
				self.x_val,
				title=self.x_title,
				bin=alt.Bin(maxbins=self.x_max, minstep=1),
				scale=alt.Scale(zero=True)
			),

			alt.Y(
				self.y_val,
				title=self.y_title
			),
			alt.Color(
				self.color_by,
				scale=alt.Scale(
					scheme=self.color_scheme
				)
			),
			tooltip=self.tooltip
		)

		return chart

	# create chart for epimap plugin
	def make_epimap(self):
		chart = alt.Chart(self.df).mark_circle(size=self.mark_size).encode(
		alt.X(
			self.x_val,
			title=self.x_title,
			axis=alt.Axis(
				values=self.x_axis_ticks
				)
			),
		alt.Y(
			self.y_val,
			title=self.y_title,
			scale=alt.Scale(type="log", reverse=True),
			axis=alt.Axis(
				values=self.y_axis_ticks
				)
			),
		alt.Color(
			self.color_by,
			scale=alt.Scale(
				scheme=self.color_scheme
				)
			),

		tooltip=self.tooltip
		)

		return chart

	# create heatmap chart for scatter zscatter plugins
	def make_heatmap_scatter(self):
		chart = alt.Chart(self.df).mark_rect().encode(
			alt.X(self.x_val, title=self.title),
			alt.X2(self.x2_val),
			alt.Y(self.y_val, title=self.y_title),
			alt.Y2(self.y2_val),
			alt.Color(self.color_by, scale = alt.Scale(
										scheme=self.color_scheme),
										legend=alt.Legend(title=self.legend_title)
										)
		)

		return chart

	# create scatterplot for scatter plugin
	def make_scatterplot(self):
		chart = alt.Chart(
			self.df
		).mark_circle(
			opacity=1, size=100, color="#009E73"
		).encode(
			alt.X(self.x_val, title = self.x_title),
			alt.Y(self.y_val, title = yTitle),
			tooltip=self.tooltip
		)

		return chart

	# create scatterplot for mutantScatters function in scatter plugin and 
	def make_mutant_scatterplot(self):
		chart = alt.Chart(self.df).mark_circle(size=self.mark_size).encode(
		x=alt.X(self.x_val, title=self.x_title),
		y=alt.Y(self.y_val, title=self.y_title),
		color=alt.Color(
			self.color_by,
			scale=alt.Scale(range=self.color_range
			)
		),
		tooltip=self.tooltip
		)

	def make__zscatterplot(self):
		alt.Chart(highlight_df).mark_point(
				filled=True, size=60
			).encode(
				x=alt.X("x:Q", title=pair[0]),
				y=alt.Y("y:Q", title=pair[1]),
				color=alt.Color(
					"highlight:N",
					scale=alt.Scale(range=[
						"#E69F00", "#56B4E9", "#009E73",
						"#F0E442", "#0072B2", "#D55E00",
						"#CC79A7"
					]),
					legend=alt.Legend(title="Significant Taxa")
				),
				# https://github.com/altair-viz/altair/issues/1181
				shape=alt.Shape(
					"highlight:N",
					legend=None
				),
				tooltip="tooltip"
			)















