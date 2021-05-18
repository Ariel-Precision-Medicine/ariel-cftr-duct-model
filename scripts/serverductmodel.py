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

def generate_source_array(input_data, key):
	'''
	Construct Arrays from Duct Model System Equations to pass along to Bokeh's ColumnDataSource

	Parameters
	----------
	input_data : dict
		Basic model parameters. New cell class with parameters instantiated each time function is called.
	key : str
		String referring to ion concentration in question. Used as a key to pull from generated nested arrays.

	Returns
	-------
	choice_array : array
		Selected array denoted by key provided of length len(time)
	'''
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
# variant_ops = pd.read_csv('cutting_variant_data.csv')['Variant']
# variant_dict = dict()
# for i in range(len(variant_ops)):
# 	variant_dict[pd.read_csv('cutting_variants.csv')['Variant'][i]] = pd.read_csv('cutting_variants.csv')['wt_func'][i]

def remove_standard_error(text_string):
	'''
	Uses unicode character code to remove standard error (+/-) icon from number

	Parameters
	----------
	text_string : str
		Example -> "132.4 +/- 4.5"

	Returns
	-------
	float
		Example -> 132.4
	'''
	return float(text_string[0:text_string.find(u'\u00b1')].strip())

variant_ops = pd.read_csv('cutting_variant_data.csv')['Variant']
variant_dict = dict()
for i in range(len(variant_ops)):
	residual = remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Residual'][i])
	ivo = remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Ivocaftor'][i])
	lum = remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Lumacaftor'][i])
	combination = remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Ivocaftor and Lumacaftor'][i])
	individual_dict = {'Residual':residual, 'Ivocaftor':ivo, 'Lumacaftor':lum, 'Combination Therapy':combination}
	variant_dict[pd.read_csv('cutting_variant_data.csv')['Variant'][i]] = individual_dict

def process_var_impact(var1, var2, therapeutic):
	'''
	Match variant impact to sweat chloride dictionary from Cutting et al

	Parameters
	----------
	var1 : str
		CFTR variant name to be used as key in Cutting data dictionary
	var2 : str
		CFTR variant name to be used as key in Cutting data dictionary
	therapeutic: str
		Therapeutic option selected on front end (if any) to examine correct column in data dictionary

	Returns
	-------
	output : float
		Average of variant impacts on chloride conductance in the presence of selected therapeutic

	'''
	vardict = {var1:None, var2:None}
	varlist = [var1, var2]
	for item in varlist:
		if item == None:
			vardict[item] = 100
		if item == 'Wild Type':
			vardict[item] = 100
		else:
			if therapeutic == 'None':
				vardict[item] = variant_dict[item]['Residual']
			if therapeutic == '10uM Ivocaftor':
				vardict[item] = variant_dict[item]['Ivocaftor']
			if therapeutic == '6uM Lumacaftor':
				vardict[item] = variant_dict[item]['Lumacaftor']
			if therapeutic == '10uM Ivocaftor + 6uM Lumacaftor':
				vardict[item] = variant_dict[item]['Combination Therapy']
	output = np.mean([vardict[var1], vardict[var2]])/100
	return output

def process_smokers(smoking_status):
	'''
	Assign Smoking Penalty to ion fluxes in data dictionary.
	smokingDictionary data estimated from Bynum et al. 1972, see supporting papers.

	Parameters
	----------
	smoking_status : str
		Description of smoking behavior (via pack-years). Five options.

	'''
	smokingDictionary = {'Non-Smoker': 1, 'Current Light Smoker': 5/11, 
						 'Current Heavy Smoker': 8/11, 'Past Light Smoker': 8.5/11,
						 'Past Heavy Smoker': 7/11}
	if smoking_status != 'Non-Smoker' and smoking_status != None:
		input_data['smoke_adj'] = smokingDictionary[smoking_status]
	else:
		input_data['smoke_adj'] = 1

def process_alcohol(alcohol_status):
	'''
	Assign Drinking Penalty to ion fluxes in data dictionary.
	Can be made more granular with the introduction of additional drinking data.

	Parameters
	----------
	alcohol_status : str
		Description of drinking behavior. Two options.

	'''
	if alcohol_status != 'Non-Drinker' and alcohol_status != None:
		input_data['alcohol_adj'] = 0.8
	else:
		input_data['alcohol_adj'] = 1

def process_therapeutics(therapeutics_status):
	'''
	Update therapeutics from front-end input

	Parameters
	----------
	therapeutics_status : str
		Current value of therapeutics radio button on front-end
	'''
	if therapeutics_status != None:
		if therapeutics_status == 'None':
			input_data['therapeutics'] = 'None'
		elif therapeutics_status == '10uM Ivocaftor':
			input_data['therapeutics'] = 'Ivocaftor'
		elif therapeutics_status == '6uM Lumacaftor':
			input_data['therapeutics'] = 'Lumacaftor'
		elif therapeutics_status == '10uM Ivocaftor + 6uM Lumacaftor':
			input_data['therapeutics'] = 'Combination Therapy'
	else:
		input_data['therapeutics'] = 'None'

def assert_WT(value):
	'''
	Catch Edge Case when WT input = None

	Parameters
	----------
	value : str
		Current string pulled from variant dropdown menu suggestion

	'''
	if value == None:
		return 'Wild Type'
	else:
		return value

# Callback Functions
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

def callback_therapeutics():
	therapeutics_radio = widgets['therapeutics']
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	index = therapeutics_radio.active
	if index != None:
		therapeutics_status = therapeutics_radio.labels[index]
		if therapeutics_status != None:
			process_therapeutics(therapeutics_status)
			input_data['variant_adj'] = process_var_impact(var1, var2, therapeutics_status)
	update_data()

def callback_luminal_antiporter():
	luminal_antiporter_radio = widgets['lum_ap']
	index =  luminal_antiporter_radio.active
	if index != None:
		lum_ap_status =  luminal_antiporter_radio.labels[index]
		if lum_ap_status == 'Off':
			input_data['ap_status'] = False
		if lum_ap_status == 'On':
			input_data['ap_status'] = True
	update_data()

def callback_basolateral_antiporter():
	luminal_antiporter_radio = widgets['baso_ap']
	index =  luminal_antiporter_radio.active
	if index != None:
		lum_ap_status =  luminal_antiporter_radio.labels[index]
		if lum_ap_status == 'Off':
			input_data['apb_status'] = False
		if lum_ap_status == 'On':
			input_data['apb_status'] = True
	update_data()

def callback_var1(attr, old, new):
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	therapeutic = input_data['therapeutics']
	input_data['variant_adj'] = process_var_impact(var1, var2, therapeutic)
	update_data()

def callback_var2(attr, old, new):
	var1, var2 = assert_WT(widgets['Variant1'].value), assert_WT(widgets['Variant2'].value)
	callback_therapeutics()
	therapeutic = input_data['therapeutics']
	input_data['variant_adj'] = process_var_impact(var1, var2, therapeutic)
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
widgets['therapeutics'].on_change('active', lambda attr, old, new: callback_therapeutics())
widgets['lum_ap'].on_change('active', lambda attr, old, new: callback_luminal_antiporter())
widgets['baso_ap'].on_change('active', lambda attr, old, new: callback_basolateral_antiporter())

# Plot Styling

# Bicarb Plot
bi_legend = Legend(items=[])
bi_plot = figure(plot_height=400, plot_width=600,  title='Patient Information')
bi_plot.add_layout(bi_legend, 'right')
bi_plot.legend.click_policy = 'hide'
bi_plot.line('time', 'pt_bi_l', source=source, line_width=3, line_color = '#34344A', legend="Patient Luminal HCO3-")
bi_plot.line('time', 'pt_bi_i', source=source, line_width=3, line_color = '#7FE0CB', legend="Patient Intra HCO3-")
bi_plot.line('time', 'wt_bi_l', source=source, line_width=3, line_color = '#34344A', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Luminal HCO3-")
bi_plot.line('time', 'wt_bi_i', source=source, line_width=3, line_color = '#7FE0CB', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Intracellular HCO3-")
bi_plot.legend.label_text_font = 'gilroy'
bi_plot.title.text_font = 'gilroy'
bi_plot.title.text_font_style = 'bold'
bi_plot.yaxis.axis_label = 'Concentration (mM)'
bi_plot.yaxis.axis_label_text_font = 'gilroy'
bi_plot.yaxis.axis_label_text_font_style = 'normal'
bi_plot.xaxis.axis_label = 'Time (min)'
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
cl_plot.yaxis.axis_label = 'Concentration (mM)'
cl_plot.yaxis.axis_label_text_font = 'gilroy'
cl_plot.yaxis.axis_label_text_font_style = 'normal'
cl_plot.xaxis.axis_label = 'Time (min)'
cl_plot.xaxis.axis_label_text_font = 'gilroy'
cl_plot.xaxis.axis_label_text_font_style = 'normal'



tabs = Tabs(tabs=[Panel(child=bi_plot, title='Bicarbonate Transport'), Panel(child=cl_plot, title='Chloride Transport')])
layout = row(input_column, tabs)

curdoc().add_root(layout)





