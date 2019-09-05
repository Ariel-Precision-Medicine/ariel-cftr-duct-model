# OOP_DUCT_MODEL.py

# For Ariel Precision Medicine
# CFTR Expert System
### Confidential ###

from bokeh.layouts import column, row, gridplot, grid
from bokeh.plotting import show
from bokeh_plotting import run_model_CFTR, graph_CFTR, graph_xd_demo, patient_plot_CFTR
from bokeh.models import Panel, Tabs
from dcw_duct_model import init_cond
from bokeh.models.widgets import Dropdown, CheckboxButtonGroup, Select, Button, Div, RadioButtonGroup, TextInput, RadioGroup
import pandas as pd

import copy

class Duct_Cell():
	# Duct cell in the pancreas
	def __init__(self):
		# Provides initial cell conditions
		self.input_dict = {'g_bi': 0.2, 'g_cl': 1, 'zeta': 0.05,
						  'kbi': 1, 'kcl': 10, 'gnbc': 2, 'gapl': 0.25,
						  'gapbl': 0.005, 'nb': 140, 'bb': 22, 'cb': 130,
						  'ni': 14, 'bi0': 15, 'buf': 0.1, 'chi': 1,
						  'bi': 15, 'ci': 60, 'bl': 32, 'gcftron': 1, 
						  'gcftrbase': 0.00007, 'ek': -0.085, 'gk': 1, 
						  'cap': 1, 'gnak': 3.125, 'np0': 25, 'epump': -0.2,
						  'ionstr': 160, 'gnaleak': 0.4, 'jac': 0.025, 
						  'rat': 0.25, 'variant_adj': None, 'vr': 0.1, 'apb_status': False,
						  'ap_status': False, 'gcftr': 0.00007, 'smoke_adj': None,
						  'alcohol_adj': None, 'therapeutics': 'None'}

		self.graphs = {'WT CFTR': None, 'Variants CFTR': None, 'Variants & Smoking CFTR': None}

		self.data = {'WT CFTR': None, 'Variants CFTR': None, 'Variants & Smoking CFTR': None}

	def report_influences(self):
		# Printout of influences added to model
		influences = [('CFTR variant','variant_adj'),
						('Smoking','smoke_adj'), ('Alcohol use', 'alcohol_adj')]
		for item in influences:
			# No changes have been implemented
			if self.input_dict[item[1]] == None:
				print(item[0] + ' factor not yet added to model')
			else:
				print(item[0] + ' factor implemented, input criteria: ' 
					+ str(self.input_dict[item[1]]))

	def add_variant_influence(self, adj):
		# Adjust the chloride conductance across membrane (% of WT fxn)
		self.input_dict['variant_adj'] = adj

	def add_smoking_influence(self, adj):
		# Adjust the functionality of CFTR (% of WT fxn)
		self.input_dict['smoke_adj'] = adj

	def add_alcohol_influence(self, adj):
		# Adjust for the influence of alcohol consumption
		self.input_dict['alcohol_adj'] = adj

	def generate_WT_graphs(self):
		self.graphs['WT CFTR'] = graph_CFTR(run_model_CFTR(init_cond, 20000, 120000, 200000),'WT_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT)')
		show(self.graphs['WT CFTR'][0])

	def generate_CFTR_graphs(self):
		input_dict = copy.deepcopy(self.input_dict)
		if self.input_dict['variant_adj'] == None:
			# WT Function
			if self.input_dict['smoke_adj'] == None and self.input_dict['alcohol_adj'] == None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'WT_plot', 'Duct Modeling Differential Equation Analysis (WT)')
			# Only Smoking
			if self.input_dict['smoke_adj'] != None and self.input_dict['alcohol_adj'] == None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Smoking_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT + Smoking)')
			# Only Alcohol
			if self.input_dict['smoke_adj'] == None and self.input_dict['alcohol_adj'] != None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT + Alcohol)')

		if self.input_dict['variant_adj'] != None:
			# Only Variant Input
			if self.input_dict['smoke_adj'] == None and self.input_dict['alcohol_adj'] == None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants)')
			# Variants and Smoking
			if self.input_dict['smoke_adj'] != None and self.input_dict['alcohol_adj'] == None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_And_Smoking_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Smoking)')
			# Variants and Alcohol
			if self.input_dict['smoke_adj'] == None and self.input_dict['alcohol_adj'] != None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_And_Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Alcohol)')
			# Variants, Smoking, and Alcohol
			if self.input_dict['smoke_adj'] != None and self.input_dict['alcohol_adj'] != None:
				self.graphs['Patient'] = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_Smoking_and_Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Smoking + Alcohol)')

	def generate_xd_demo_data(self):
		# Run Mmdel with no influences
		self.data['WT CFTR'] = run_model_CFTR(init_cond, 20000, 120000, 200000)
		# Run model with variant influences
		if self.input_dict['variant_adj'] != None:
			variant_input_dict = copy.deepcopy(self.input_dict)
			variant_input_dict['smoke_adj'] = None
			variant_input_dict['alcohol_adj'] = None
			self.data['Variants CFTR'] = run_model_CFTR(variant_input_dict, 20000, 120000, 200000)
		# Run model with smoking and variants
		if self.input_dict['variant_adj'] != None and self.input_dict['smoke_adj'] != None:
			variant_and_smoking_input_dict = copy.deepcopy(self.input_dict)
			variant_and_smoking_input_dict['alcohol_adj'] = None
			self.data['Variants & Smoking CFTR'] = run_model_CFTR(variant_and_smoking_input_dict, 20000, 120000, 200000)

	def generate_xd_graphs(self):
		self.generate_xd_demo_data()
		# Call graphing function
		graph_xd_demo(self.data['WT CFTR'], self.data['Variants CFTR'], self.data['Variants & Smoking CFTR'])

	def generate_tabbed_output(self):
		self.patient_graph()
		tab_bicarb, tab_chloride = self.graphs['Patient']
		tab_bicarb.plot_width = 800
		tab_chloride.plot_width = 800
		panel_bicarb = Panel(child=tab_bicarb, title="Bicarbonate Transport")
		panel_chloride = Panel(child=tab_chloride, title="Chloride Transport")
		tabs = Tabs(tabs=[ panel_bicarb, panel_chloride ], width=1000)
		var_menu = self.gen_var_menu()
		widget_column, widgets = self.process_widgets(var_menu)
		layout = grid([[widget_column,tabs]])
		return layout, widgets, tab_bicarb, tab_chloride

	def gen_var_menu(self):
		self.variant_ops = pd.read_csv('cutting_variant_data.csv')['Variant']
		self.variant_dict = dict()
		for i in range(len(self.variant_ops)):
			residual = self.remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Residual'][i])
			ivo = self.remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Ivocaftor'][i])
			lum = self.remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Lumacaftor'][i])
			combination = self.remove_standard_error(pd.read_csv('cutting_variant_data.csv')['Ivocaftor and Lumacaftor'][i])
			individual_dict = {'Residual':residual, 'Ivocaftor':ivo, 'Lumacaftor':lum, 'Combination Therapy':combination}
			self.variant_dict[pd.read_csv('cutting_variant_data.csv')['Variant'][i]] = individual_dict
		var_menu = []
		for key in self.variant_dict:
			var_menu.append(key)
		var_menu.insert(0, 'Wild Type')
		return var_menu

	def remove_standard_error(self, text_string):
		return text_string[0:text_string.find(u'\u00b1')].strip()

	def process_widgets(self, var_menu):
		select_cftr1 = Select(title="CFTR Variant 1", value=None, options=var_menu, width = 100)
		select_cftr2 = Select(title="CFTR Variant 2", value=None, options=var_menu, width = 100)
		select_row = row(select_cftr1, select_cftr2, width = 500)
		text_input = TextInput(width=200)
		smoking_radio = RadioGroup(labels=["Non-Smoker", "Current Light Smoker", "Current Heavy Smoker", "Past Light Smoker", "Past Heavy Smoker"], active = 0)
		drinking_radio = RadioButtonGroup(labels=["Non-Drinker", "Drinker"], active = 0)
		therapeutic_radio = RadioButtonGroup(labels=["None", "10uM Ivocaftor", "6uM Lumacaftor", "10uM Ivocaftor + 6uM Lumacaftor"], active = 0)
		luminal_antiporter_radio = RadioButtonGroup(labels=["Off", "On"], active = 0)
		basolateral_antiporter_radio = RadioButtonGroup(labels=["Off", "On"], active = 0)
		pt_div = Div(text = """<b>Patient Identifier</b>""", width=110, height=25, style={'font-family':'gilroy'})
		ariel_div = Div(text = """<b>Duct Model </b> [Property of <a 
						href="https://www.arielmedicine.com/" target="_blank">Ariel Precision 
						Medicine</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		smoking_div = Div(text = """<b>Smoking Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1412225/"
						target="_blank">PMID 5036091</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		drinking_div = Div(text = """<b>Drinking Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4353632/"
						target="_blank">PMID 25447846</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		variants_div = Div(text = """<b>Patient CFTR Variants </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124440/"
						target="_blank">PMID 30046002</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		therapeutics_div = Div(text = """<b>Therapeutics </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124440/"
						target="_blank">PMID 30046002</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		lum_ap_div = Div(text = """<b>Luminal Antiporter Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pubmed/15257112"
						target="_blank">PMID 15257112</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		baso_ap_div = Div(text = """<b>Basolateral Antiporter Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pubmed/15257112"
						target="_blank">PMID 15257112</a>]""", width=500, height=25, style={'font-family':'gilroy'})
		#image_div = Div(width = 200, height = 100, text="<img src='https://pancreapedia.org/sites/default/files/Figure1%20copy.jpg'>")
		save_button = Button(label='Save Patient Data', button_type='warning')
		input_row = row(pt_div, text_input, width = 500)
		widget_column = column(ariel_div, variants_div, select_row,
								therapeutics_div, therapeutic_radio,
								smoking_div, smoking_radio,
								drinking_div, drinking_radio,
								lum_ap_div, luminal_antiporter_radio,
								baso_ap_div, basolateral_antiporter_radio)
								#image_div
		widgets = {'Variant1': select_cftr1, 'Variant2': select_cftr2,
					'smoking_status': smoking_radio, 'alcohol_status': drinking_radio,
					'pt_id':text_input, 'therapeutics':therapeutic_radio, 'lum_ap':luminal_antiporter_radio,
					'baso_ap':basolateral_antiporter_radio}
		# widgets = [select_cftr1, select_cftr2, smoking_radio, drinking_radio]
		# for item in widgets:
		# 	item.on_change('value', self.user_inputs)
		return widget_column, widgets

	# def user_inputs(self):
	# 	self.variant_adj = process_variants(self.widgets['Variant1'].value, self.widgets['Variant2'].value)
	# 	self.smoke_adj = process_smoking(self.widgets['smoking_status'].value)
	# 	self.alcohol_adj = process_alcohol(self.widgets['alcohol_status'].value)

	# def process_variants(var1, var2):
	# 	print(var1, var2)

	# def process_smoking(smoking_status):
	# 	print(smoking_status)

	# def process_alcohol(alcohol_status):
	# 	print(alcohol_status)

	def patient_graph(self):
		input_dict = copy.deepcopy(self.input_dict)
		p1, p2 = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_Smoking_and_Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Smoking + Alcohol)')
		wt_results = run_model_CFTR(init_cond, 20000, 120000, 200000)
		self.graphs['Patient'] = patient_plot_CFTR(p1, p2, wt_results, None, None)

	def patient_graph_server(input_dict):
		p1, p2 = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_Smoking_and_Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Smoking + Alcohol)')
		wt_results = run_model_CFTR(init_cond, 20000, 120000, 200000)
		return patient_plot_CFTR(p1, p2, wt_results, None, None)


#cell1 = Duct_Cell()
#cell1.generate_CFTR_graphs()
#show(cell1.graphs['Patient'][1])

