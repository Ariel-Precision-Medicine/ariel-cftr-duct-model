# OOP_DUCT_MODEL.py

# For Ariel Precision Medicine
# CFTR Expert System
### Confidential ###

from bokeh.layouts import column, row, gridplot, grid
from bokeh.plotting import show
from bokeh_plotting import run_model_CFTR, graph_CFTR, graph_xd_demo
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

	def generate_CFTR_graphs(self):
		rows = []
		self.graphs['WT CFTR'] = graph_CFTR(run_model_CFTR(init_cond, 20000, 120000, 200000),'WT_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT)')
		rows.append(self.graphs['WT CFTR'])

		if self.input_dict['variant_adj'] != None:
			variant_input_dict = copy.deepcopy(self.input_dict)
			variant_input_dict['smoke_adj'] = None
			variant_input_dict['alcohol_adj'] = None
			self.graphs['Variants CFTR'] = graph_CFTR(run_model_CFTR(variant_input_dict, 20000, 120000, 200000), 'Variants_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT + Variants)')
			rows.append(self.graphs['Variants CFTR'])

		if self.input_dict['variant_adj'] != None and self.input_dict['smoke_adj'] != None:
			variant_and_smoking_input_dict = copy.deepcopy(self.input_dict)
			variant_and_smoking_input_dict['alcohol_adj'] = None
			self.graphs['Variants & Smoking CFTR'] = graph_CFTR(run_model_CFTR(variant_and_smoking_input_dict, 20000, 120000, 200000), 'Variants_And_Smoking_CFTR_plot', 'Duct Modeling Differential Equation Analysis (WT + Variants + Smoking)')
			rows.append(self.graphs['Variants & Smoking CFTR'])

		# make a grid
		grid = gridplot(rows, plot_width = 550, plot_height = 350, sizing_mode = 'scale_width')
		
		# show the results
		show(grid)

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

		self.variant_ops = pd.read_csv('cutting_variants.csv')['variant']
		self.variant_dict = dict()
		for i in range(len(self.variant_ops)):
			self.variant_dict[pd.read_csv('cutting_variants.csv')['variant'][i]] = pd.read_csv('cutting_variants.csv')['wt_func'][i]
		menu = []
		for key in self.variant_dict:
			menu.append(key)

		# dropdown = Dropdown(label="CFTR Variant List", button_type="danger", menu=menu)
		menu.insert(0, '-')
		select = Select(title="CFTR Variants", value=None, options=menu, width = 100)
		button_var_add = Button(label="Add Variant")
		button_var_remove = Button(label="Remove Variant")
		current_variants_string = """see variants here,see variants here,see variants here,see variants here,see variants here"""
		div = Div(text=current_variants_string,width=300, height=50)
		select_row = row(select, div, width = 400)
		var_row = row(button_var_add, button_var_remove, width = 400)
		smoking_radio = RadioButtonGroup(labels=["Non-Smoker", "Light Smoker", "Heavy Smoker", "Past Smoker"])
		drinking_radio = RadioButtonGroup(labels=["Non-Drinker", "Light Drinker", "Heavy Drinker"])
		ariel_div = Div(text = """<b>Duct Model </b> [Property of <a 
						href="https://www.arielmedicine.com/" target="_blank">Ariel Precision 
						Medicine</a>]""", width=400, height=25)
		smoking_div = Div(text = """<b>Patient Smoking Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1412225/"
						target="_blank">PMID 5036091</a>]""", width=400, height=25)
		drinking_div = Div(text = """<b>Patient Drinking Status </b> [Source: <a 
						href="https://www.google.com/"
						target="_blank">PMID Needs Updated</a>]""", width=400, height=25)
		variants_div = Div(text = """<b>Patient CFTR Status </b> [Source: <a 
						href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124440/"
						target="_blank">PMID 30046002</a>]""", width=400, height=25)
		input_column = column(ariel_div, variants_div, select_row,
								var_row, smoking_div, smoking_radio,
								drinking_div, drinking_radio)
		l = grid([[input_column,tabs]])
		show(l)




wt_cell = Duct_Cell()
wt_cell.report_influences()
# wt_cell.add_smoking_influence(0.50)
wt_cell.add_variant_influence(0.266)
wt_cell.report_influences()
# wt_cell.generate_CFTR_graphs()
# wt_cell.generate_xd_graphs()
wt_cell.generate_tabbed_output()
 




