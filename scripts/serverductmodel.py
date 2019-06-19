from oop_duct_model import Duct_Cell
from bokeh.io import curdoc, export_png
from bokeh.plotting import show, figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh_plotting import run_model_CFTR
from bokeh.models import ColumnDataSource, Legend, Tabs, Panel
from bokeh.models.widgets import Dropdown, CheckboxButtonGroup, Select, Button, Div, RadioButtonGroup, TextInput
import pandas as pd
import numpy as np
from dcw_duct_model import init_cond
from bokeh.layouts import row
import copy

# Script to generate duct model instance via server
# Callbacks connect to Python logic and return
# new information to be graphed

# Instantiate New Cell for Patient
wt_cell = Duct_Cell()
input_data = wt_cell.input_dict

init = copy.deepcopy(init_cond)

# Construct Arrays from Duct Model System Equations
def generate_source_array(input_data, key):
	choice_array = run_model_CFTR(input_data, 20000, 120000, 200000)[0][key]
	if key == 'time':
		choice_array /=  20000
	return choice_array

# Individual graphing arrays for ColumnDataSource
wt_bi_l = generate_source_array(init, 'bl')
wt_bi_i = generate_source_array(init, 'bi')
pt_bi_l = generate_source_array(input_data, 'bl')
pt_bi_i = generate_source_array(input_data, 'bi')
wt_cl_l = 160 - wt_bi_l
wt_cl_i = generate_source_array(init, 'ci')
pt_cl_l = 160 - pt_bi_l
pt_cl_i = generate_source_array(input_data, 'ci')
time = generate_source_array(input_data, 'time')

# Construct CDS to store data in server
source = ColumnDataSource(data=dict(wt_bi_l = wt_bi_l,
									wt_bi_i = wt_bi_i,
									pt_bi_l = pt_bi_l,
									pt_bi_i = pt_bi_i,
									wt_cl_l = wt_cl_l,
									wt_cl_i = wt_cl_i,
									pt_cl_l = pt_cl_l,
									pt_cl_i = pt_cl_i,
									time = time))

# Build User Input Column Widgets
input_column, widgets = wt_cell.process_widgets(wt_cell.gen_var_menu())

# Build Variant Resource Dictionary
variant_ops = pd.read_csv('cutting_variants.csv')['variant']
variant_dict = dict()
for i in range(len(variant_ops)):
	variant_dict[pd.read_csv('cutting_variants.csv')['variant'][i]] = pd.read_csv('cutting_variants.csv')['wt_func'][i]

# Match variant impact to sweat chloride dictionary
def process_var_impact(var1, var2):
	if var1 == 'Wild Type':
		var1_val = 100
	else:
		var1_val = variant_dict[var1]
	if var2 == 'Wild Type':
		var2_val = 100
	else:
		var2_val = variant_dict[var2]
	return np.mean([var1_val, var2_val])/100

# Assign Smoking Penalty
def process_smokers(smoking_status):
	if smoking_status != 'Non-Smoker' and smoking_status != None:
		input_data['smoke_adj'] = 0.75
	else:
		input_data['smoke_adj'] = 1

# Assign Drinking Penalty
def process_alcohol(alcohol_status):
	if alcohol_status != 'Non-Drinker' and alcohol_status != None:
		input_data['alcohol_adj'] = 0.8
	else:
		input_data['alcohol_adj'] = 1

# Catch Edge Case when WT input == None
def assert_WT(value):
	if value == None:
		return 'Wild Type'
	else:
		return value

# Callback Functions
def callback_var1(attr, old, new):
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	input_data['variant_adj'] = process_var_impact(var1, var2)
	update_data()

def callback_var2(attr, old, new):
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	input_data['variant_adj'] = process_var_impact(var1, var2)
	update_data()

def callback_smoking():
	smoking_radio = widgets['smoking_status']
	index = smoking_radio.active
	smoking_status = smoking_radio.labels[index]
	process_smokers(smoking_status)
	update_data()

def callback_drinking():
	drinking_radio = widgets['alcohol_status']
	index = drinking_radio.active
	alcohol_status = drinking_radio.labels[index]
	process_alcohol(alcohol_status)
	update_data()

def update_data():
	pt_bi_l = generate_source_array(input_data, 'bl')
	pt_bi_i = generate_source_array(input_data, 'bi')
	wt_bi_l = generate_source_array(init, 'bl')
	wt_bi_i = generate_source_array(init, 'bi')
	wt_cl_l = (160 - wt_bi_l)
	wt_cl_i = generate_source_array(init, 'ci')
	pt_cl_l = 160 - pt_bi_l
	pt_cl_i = generate_source_array(input_data, 'ci')
	time = generate_source_array(input_data, 'time')

	source.data = dict( wt_bi_l = wt_bi_l,
						wt_bi_i = wt_bi_i,
						pt_bi_l = pt_bi_l,
						pt_bi_i = pt_bi_i,
						wt_cl_l = wt_cl_l,
						wt_cl_i = wt_cl_i,
						pt_cl_l = pt_cl_l,
						pt_cl_i = pt_cl_i,
						time = time)

# Event Calls
widgets['Variant1'].on_change('value', callback_var1)
widgets['Variant2'].on_change('value', callback_var2)
widgets['smoking_status'].on_change('active', lambda attr, old, new: callback_smoking())
widgets['alcohol_status'].on_change('active', lambda attr, old, new: callback_drinking())


# Plot Styling

# Bicarb Plot
bi_legend = Legend(items=[])
bi_plot = figure(plot_height=400, plot_width=600)
bi_plot.add_layout(bi_legend, 'right')
bi_plot.legend.click_policy = 'hide'
bi_plot.line('time', 'pt_bi_l', source=source, line_width=3, line_color = '#34344A', legend="Patient Luminal HCO3-")
bi_plot.line('time', 'pt_bi_i', source=source, line_width=3, line_color = '#7FE0CB', legend="Patient Intra HCO3-")
bi_plot.line('time', 'wt_bi_l', source=source, line_width=3, line_color = '#34344A', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Luminal HCO3-")
bi_plot.line('time', 'wt_bi_i', source=source, line_width=3, line_color = '#7FE0CB', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Intracellular HCO3-")
bi_plot.legend.label_text_font = 'gilroy'
bi_plot.title.text_font = 'gilroy'
bi_plot.title.text_font_style = 'bold'
bi_plot.yaxis.axis_label_text_font = 'gilroy'
bi_plot.yaxis.axis_label_text_font_style = 'normal'
bi_plot.xaxis.axis_label_text_font = 'gilroy'
bi_plot.xaxis.axis_label_text_font_style = 'normal'

# Chloride Plot
cl_legend = Legend(items=[])
cl_plot = figure(plot_height=400, plot_width=600,  title='Patient Information')
cl_plot.add_layout(cl_legend, 'right')
cl_plot.legend.click_policy = 'hide'
cl_plot.line('time', 'pt_cl_l', source=source, line_width=3, line_color = '#34344A', legend="Patient Luminal Cl-")
cl_plot.line('time', 'pt_cl_i', source=source, line_width=3, line_color = '#7FE0CB', legend="Patient Intra Cl-")
cl_plot.line('time', 'wt_cl_l', source=source, line_width=3, line_color = '#34344A', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Luminal Cl-")
cl_plot.line('time', 'wt_cl_i', source=source, line_width=3, line_color = '#7FE0CB', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Intracellular Cl-")
cl_plot.legend.label_text_font = 'gilroy'
cl_plot.legend.background_fill_color = '#F4F1E1'
cl_plot.legend.background_fill_alpha = 0.25
cl_plot.title.text_font = 'gilroy'
cl_plot.title.text_font_style = 'bold'
cl_plot.yaxis.axis_label_text_font = 'gilroy'
cl_plot.yaxis.axis_label_text_font_style = 'normal'
cl_plot.xaxis.axis_label_text_font = 'gilroy'
cl_plot.xaxis.axis_label_text_font_style = 'normal'

tabs = Tabs(tabs=[Panel(child=bi_plot, title='Bicarbonate Transport'), Panel(child=cl_plot, title='Chloride Transport')])
layout = row(input_column, tabs)

curdoc().add_root(layout)





