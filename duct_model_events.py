# duct_model_events.py

'''
Callback events linking Pancreatic Duct Model GUI to logic in
DCW duct model script

Use this as a tool to observe mechanistic view of duct cell

Developed by Ariel Precision Medicine
'''

import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5 import QtCore, QtGui, QtWidgets
import pandas as pd
import numpy as np

# Duct Model MainWindow
from duct_model_gui import Ui_MainWindow

# Import mathematical modeling functions
from DCW_duct_model import *
from DCW_duct_graphing_functions import *
from math import floor

class MainWindow(QMainWindow, Ui_MainWindow):

	def __init__(self):
		print('In MainWindow init')
		QMainWindow.__init__(self)
		self.setupUi(self)

		# Dictionary with Default Duct Model Parameters
		self.duct_model_default_inputs = init_cond

		# Graphing Options
		self.graph_ops = ['-', 'GCFTR', 'Volume Ratios', 'Antiporters']

		# Variant Options
		self.variant_ops = pd.read_csv('cutting_variants.csv')['variant']
		self.variant_dict = dict()
		for i in range(len(self.variant_ops)):
			self.variant_dict[pd.read_csv('cutting_variants.csv')['variant'][i]] = pd.read_csv('cutting_variants.csv')['wt_func'][i]
		self.comboBox_variant_input_area.addItem('-')
		self.comboBox_variant_input_area.addItems(self.variant_ops)

		self.wt_total = 100

		# List of Patient Variants
		self.variant_set = set()
		self.patient_variants = []

		# Connect User Actions
		self.connect_elements()

	def string_of_variants(self, variant_list):
		output_string = ''
		for item in variant_list:
			output_string += (item + ' %f' % round(self.variant_dict[item],2) + '%' + '\n')
		return output_string

	def connect_elements(self):
		'''
		Connect GUI widgets to their interaction and update functions
		'''
		self.pushButton_add_variant.clicked.connect(self.add_variant)
		self.pushButton_remove_variant.clicked.connect(self.remove_variant)

	def update_wt_func(self):
		self.progressBar_wt_functionality.setValue(self.wt_total)

	def collect_variant_impact(self):
		impact_percentage = 0
		variant_impacts = []

		for item in self.patient_variants:
			variant_impacts.append(float(self.variant_dict[item]))
		impact_percentage = np.mean(variant_impacts)
		print(impact_percentage)
		if impact_percentage > 100:
			impact_percentage = 100
			print('Variants perform above expected WT functionality by %f percent' % (impact_percentage-100))
		return impact_percentage

	def add_variant(self):
		variant = self.comboBox_variant_input_area.currentText()
		if variant != '-':
			if variant not in self.variant_set:
				self.patient_variants.append(variant)
				text = self.string_of_variants(self.patient_variants)
				self.plainTextEdit_patient_variants.setPlainText(text)
				self.variant_set.add(variant)
				self.wt_total = self.collect_variant_impact()
				self.update_wt_func()

	def remove_variant(self):
		variant = self.comboBox_variant_input_area.currentText()
		if variant != '-':
			if variant in self.variant_set:
				self.patient_variants.remove(variant)
				text = self.string_of_variants(self.patient_variants)
				self.plainTextEdit_patient_variants.setPlainText(text)
				self.variant_set.remove(variant)
				self.wt_total = self.collect_variant_impact()
				self.update_wt_func()

	
	def update_patient_variants(self):
		self.plainTextEdit_patient_variants.setPlainText(str(self.patient_variants))




if __name__ == '__main__':

	app = QApplication(sys.argv)
	# Duct Model GUI Window
	duct_model_app = MainWindow()
	duct_model_app.show()

	sys.exit(app.exec_())
