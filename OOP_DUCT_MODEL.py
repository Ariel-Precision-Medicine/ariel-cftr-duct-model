# OOP_DUCT_MODEL.py

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


cell1 = Duct_Cell()
cell1.report_influences()
print('\n')
cell1.add_variant_influence(0.85)
cell1.add_alcohol_influence(0.70)
cell1.add_smoking_influence(0.82)
cell1.report_influences()