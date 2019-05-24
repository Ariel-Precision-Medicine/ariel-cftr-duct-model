'''
This script tests the duct model for the Celeste Repository / Expert System.

Test functions are written for the differential equations of the duct model,
graphing logic, and variant logic.

Developed by Ariel Precision Medicine. 
'''
import unittest, math
from DCW_duct_model import *
from DCW_duct_graphing_functions import *
# TODO: Need to add test cases for GUI logic and UI

init_cond = {'g_bi': 0.2, 'g_cl': 1, 'zeta': 0.05,
			  'kbi': 1, 'kcl': 10, 'gnbc': 2, 'gapl': 0.25,
			  'gapbl': 0.005, 'nb': 140, 'bb': 22, 'cb': 130,
			  'ni': 14, 'bi0': 15, 'buf': 0.1, 'chi': 1,
			  'bi': 15, 'ci': 60, 'bl': 32, 'gcftron': 1, 
			  'gcftrbase': 0.00007, 'ek': -0.085, 'gk': 1, 
			  'cap': 1, 'gnak': 3.125, 'np0': 25, 'epump': -0.2,
			  'ionstr': 160, 'gnaleak': 0.4, 'jac': 0.025, 
			  'rat': 0.25, 'cond_adj': 1, 'vr': 0.1, 'apb_status': True,
			  'ap_status': True, 'gcftr': 0.00007}


class TestDuctModelLogic(unittest.TestCase):
	''' 
		Load Initial Conditions from 
		"A Mathematical Model of the Pancreatic Duct Cell
		Generating High Bicarbonate Concentrations in
		Pancreatic Juice" 

		https://pdfs.semanticscholar.org/b312/3f29dbb27090bac662e0fc6bb452ed450a79.pdf
	'''
	def setUp(self, *initial_conditions, **kwargs):
		for dictionary in initial_conditions:
			for key in dictionary:
				setattr(self, key, dictionary[key])

	def test_antiporter_math(self):
		# Test Fxn Math

		# Numerator = 0
		self.assertEqual(antiporter(1,1,1,1,1,1), 0)
		self.assertEqual(antiporter(1,1,2,2,1,1), 0)
		self.assertEqual(antiporter(2,1,2,1,1,1), 0)
		self.assertEqual(antiporter(1,2,1,2,1,1), 0)
		# Non-zero Numerator
		self.assertEqual(antiporter(2,1,1,2,1,1), 3/24)

	def test_eff_perm(self):
		pass

	def test_nernst_potential(self):
		# Test Fxn Math

		self.assertEqual(nernst_potential(math.e, math.e), 0)
		self.assertEqual(nernst_potential(math.e**3, math.e), 8.31451*310/96485*2)

	def test_duct_model_system(self):
		pass

class TestDuctModelGraphingFunctions(unittest.TestCase):
	''' 
		Load Initial Conditions from 
		"A Mathematical Model of the Pancreatic Duct Cell
		Generating High Bicarbonate Concentrations in
		Pancreatic Juice" 

		https://pdfs.semanticscholar.org/b312/3f29dbb27090bac662e0fc6bb452ed450a79.pdf
	'''
	def setUp(self, *initial_conditions, **kwargs):
		for dictionary in initial_conditions:
			for key in dictionary:
				setattr(self, key, dictionary[key])

	def test_graph_generation(self):
		pass

	def test_cftr_calc_HCO3_Cl(self):
		pass

	def test_vol_rat_calc(self):
		pass

	def test_antiporters_calc(self):
		pass

	def test_calc_variant_impact(self):
		pass
		

def suite():
	suite = unittest.TestSuite()
	suite.addTest(unittest.makeSuite(TestDuctModelLogic))
	suite.addTest(unittest.makeSuite(TestDuctModelGraphingFunctions))
	return suite

if __name__ == '__main__':
	runner = unittest.TextTestRunner()
	runner.run(suite())


