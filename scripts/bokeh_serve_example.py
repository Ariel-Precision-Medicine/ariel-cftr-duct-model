from oop_duct_model import Duct_Cell
from bokeh.io import curdoc
from bokeh.plotting import show
from bokeh_plotting import run_model_CFTR
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Dropdown, CheckboxButtonGroup, Select, Button, Div, RadioButtonGroup
import pandas as pd
import numpy as np
from dcw_duct_model import init_cond


# Instantiate New Cell for Patient
wt_cell = Duct_Cell()
input_data = wt_cell.input_dict
print(input_data)

# Construct Arrays from Duct Model System Equations
def generate_source_array(input_data, key):
	choice_array = run_model_CFTR(input_data, 20000, 120000, 200000)[0][key]
	if key == 't':
		choice_array / 20000
	return choice_array

# Individual graphing arrays for ColumnDataSource
wt_bi_l = generate_source_array(init_cond, 'bl')
wt_bi_i = generate_source_array(init_cond, 'bi')
pt_bi_l = generate_source_array(input_data, 'bl')
pt_bi_i = generate_source_array(input_data, 'bi')
wt_cl_l = 160 - wt_bi_l
wt_cl_i = generate_source_array(init_cond, 'ci')
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
	return np.mean([var1_val, var2_val])

# Assign Smoking Penalty
def process_smokers(smoking_status):
	if smoking_status != 'Non-Smoker':
		input_data['smoke_adj'] = 0.75
	else:
		input_data['smoke_adj'] = 1

# Assign Drinking Penalty
def process_alcohol(alcohol_status):
	if alcohol_status != 'Non-Drinker':
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

def callback_var2(attr, old, new):
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	input_data['variant_adj'] = process_var_impact(var1, var2)

def callback_smoking():
	smoking_radio = widgets['smoking_status']
	index = smoking_radio.active
	smoking_status = smoking_radio.labels[index]
	process_smokers(smoking_status)
	print(input_data)
	print(source.data)

def callback_drinking():
	drinking_radio = widgets['alcohol_status']
	index = drinking_radio.active
	alcohol_status = drinking_radio.labels[index]
	process_alcohol(alcohol_status)
	print(input_data)
	print(source.data)

# Event Calls
widgets['Variant1'].on_change('value', callback_var1)
widgets['Variant2'].on_change('value', callback_var2)
widgets['smoking_status'].on_change('active', lambda attr, old, new: callback_smoking())
widgets['alcohol_status'].on_change('active', lambda attr, old, new: callback_drinking())

curdoc().add_root(input_column)
