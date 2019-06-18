# OOP_DUCT_MODEL.py

# For Ariel Precision Medicine
# CFTR Expert System
### Confidential ###

from bokeh.layouts import column, row, gridplot, grid
from bokeh.plotting import show
from bokeh_plotting import run_model_CFTR, graph_CFTR, graph_xd_demo, patient_plot_CFTR
from bokeh.models import Panel, Tabs
from dcw_duct_model import init_cond
from bokeh.models.widgets import Dropdown, CheckboxButtonGroup, Select, Button, Div, RadioButtonGroup
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
						  'alcohol_adj': None}

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
		self.graphs['WT CFTR'] = graph_CFTR(run_model_CFTR(init_cond, 20000, 120000, 200000),'WT_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT)')
		tab_bicarb, tab_chloride = self.graphs['WT CFTR']
		tab_bicarb.plot_width = 800
		tab_chloride.plot_width = 800
		tab_bicarb = Panel(child=tab_bicarb, title="Bicarbonate Transport")
		tab_chloride = Panel(child=tab_chloride, title="Chloride Transport")
		tabs = Tabs(tabs=[ tab_bicarb, tab_chloride ], width=1000)
		var_menu = self.gen_var_menu()
		display = self.process_widgets(tabs, var_menu)
		show(display)

	def gen_var_menu(self):
		self.variant_ops = pd.read_csv('cutting_variants.csv')['variant']
		self.variant_dict = dict()
		for i in range(len(self.variant_ops)):
			self.variant_dict[pd.read_csv('cutting_variants.csv')['variant'][i]] = pd.read_csv('cutting_variants.csv')['wt_func'][i]
		var_menu = []
		for key in self.variant_dict:
			var_menu.append(key)
		var_menu.insert(0, 'Wild Type')
		return var_menu

	def process_widgets(self, tabs, var_menu):
		select_cftr1 = Select(title="CFTR Variant 1", value=None, options=var_menu, width = 100)
		select_cftr2 = Select(title="CFTR Variant 2", value=None, options=var_menu, width = 100)
		button_var_reset = Button(label="Reset Variants", width=75)
		current_variants_string = """see variants here,see variants here,see variants here,see variants here,see variants here"""
		select_row = row(select_cftr1, select_cftr2, button_var_reset, width = 400)
		smoking_radio = RadioButtonGroup(labels=["Non-Smoker", "Light Smoker", "Heavy Smoker", "Past Smoker"])
		drinking_radio = RadioButtonGroup(labels=["Non-Drinker", "Drinker"])
		ariel_div = Div(text = """<b>Duct Model </b> [Property of <a 
						href="https://www.arielmedicine.com/" target="_blank">Ariel Precision 
						Medicine</a>]""", width=400, height=25)
		smoking_div = Div(text = """<b>Patient Smoking Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1412225/"
						target="_blank">PMID 5036091</a>]""", width=400, height=25)
		drinking_div = Div(text = """<b>Patient Drinking Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4353632/"
						target="_blank">PMID 25447846</a>]""", width=400, height=25)
		variants_div = Div(text = """<b>Patient CFTR Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124440/"
						target="_blank">PMID 30046002</a>]""", width=400, height=25)
		input_column = column(ariel_div, variants_div, select_row,
								smoking_div, smoking_radio,
								drinking_div, drinking_radio)
		layout = grid([[input_column,tabs]])
		return layout


	def fancy_graph(self):
		input_dict = copy.deepcopy(self.input_dict)
		p1, p2 = graph_CFTR(run_model_CFTR(input_dict, 20000, 120000, 200000), 'Variants_Smoking_and_Alcohol_CFTR_plot', 'Duct Modeling Differential Equation Analysis (Variants + Smoking + Alcohol)')
		wt_results = run_model_CFTR(init_cond, 20000, 120000, 200000)
		p3, p4 = patient_plot_CFTR(p1, p2, wt_results, None, None)
		show (p3)
		show (p4)

wt_cell = Duct_Cell()
wt_cell.report_influences()
wt_cell.add_smoking_influence(0.30)
wt_cell.add_variant_influence(0.266)
wt_cell.add_alcohol_influence(0.75)
wt_cell.report_influences()
wt_cell.fancy_graph()
# wt_cell.generate_WT_graphs()
# wt_cell.generate_CFTR_graphs()
# wt_cell.generate_xd_graphs()
# wt_cell.generate_tabbed_output()
 




