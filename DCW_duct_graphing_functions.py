#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 11:18:16 2018

@author: crbreze
for Ariel Precision Medicine

Translation of the Basic Pancreatic Duct Cell Modeling Script described in the
appendix of Whitcomb & Ermentrout's 'A Mathematical Model of the Pancreatic 
Duct Cell Generating High Bicarbonate Concentrations in Pancreatic Juice'
"""
from math import log
from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rnd
import copy
from DCW_duct_model import duct_model_system, init_cond

def graph_generation(graph_type, input_dict):
	filename = None
	if graph_type == 'GCFTR':
		filename = graph_CFTR(cftr_calc_HCO3_Cl(input_dict, 20000, 120000, 200000),'GCFTR')
	elif graph_type == 'Volume Ratios':
		volumes_list = [10, 5, 1, 0.5, 0.25, 0.1, 0.01]
		filename = graph_volume_ratios(vol_rat_calc(input_dict, 20000, 120000, 200000, volumes_list),'VolumeRatio')
	elif graph_type == 'Antiporters':
		filename = graph_antiporters(antiporters_calc(input_dict, 20000, 120000, 200000), 'antiporter_onoff')

	return filename

def cftr_calc_HCO3_Cl(input_dict, t_on, t_off, t_end):
	# Wrapper function to pass dictionary parameters through solve_ivp
	def wrapper_fxn(t, y):
		return duct_model_system(t, y, cond)
	# Copy Input Dictionary so Initial Dictionary is Unaltered
	cond = copy.deepcopy(input_dict)
	# Period Before gcftr opens
	cond['gcftr'] = init_cond['gcftrbase']
	y0_0 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
	state0 = solve_ivp(wrapper_fxn, [0,t_on], y0_0)
	# Period where gcftr is open
	# Turn gcftr 'on' and re-initialize boundary conditions
	cond['gcftr'] = cond['gcftron']
	cond['bi'] = state0['y'][0][-1]
	cond['bl'] = state0['y'][1][-1]
	cond['ci'] = state0['y'][2][-1]
	cond['ni'] = state0['y'][3][-1]
	y0_1 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
	state1 = solve_ivp(wrapper_fxn, [t_on,t_off], y0_1)
	# Period where gcftr closes to end
	# Turn gcftr 'on' and re-initialize boundary conditions
	cond['gcftr'] = cond['gcftrbase']
	cond['bi'] = state1['y'][0][-1]
	cond['bl'] = state1['y'][1][-1]
	cond['ci'] = state1['y'][2][-1]
	cond['ni'] = state1['y'][3][-1]
	y0_2 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
	state2 = solve_ivp(wrapper_fxn, [t_off,t_end], y0_2)

	t_graph = np.concatenate([state0['t'], state1['t'], state2['t']])
	bi_graph = np.concatenate([state0['y'][0], state1['y'][0], state2['y'][0]])
	bl_graph = np.concatenate([state0['y'][1], state1['y'][1], state2['y'][1]])
	ci_graph = np.concatenate([state0['y'][2], state1['y'][2], state2['y'][2]])
	ni_graph = np.concatenate([state0['y'][3], state1['y'][3], state2['y'][3]])

	return [[t_graph, bi_graph, bl_graph, ci_graph, ni_graph], t_on, t_off, t_end]

def graph_CFTR(model_results, filename):
	# Unpack variables
	t = model_results[0][0]
	bi = model_results[0][1]
	bl = model_results[0][2]
	ci = model_results[0][3]
	ni = model_results[0][4]
	t_on = model_results[1]
	t_off = model_results[2]
	t_end = model_results[3]

	# Upper Plot
	plt.subplot(2, 1, 1)
	plt.plot(t,bi, 'r-', label = 'HCO3-(intra)')
	plt.plot(t,bl, 'g-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.title('Duct Modeling Dif. Eq. \n GCFTR OPEN in RED')
	plt.ylabel('Bicarb Conc. (mM)')

	# Lower Plot
	plt.subplot(2, 1, 2)
	plt.plot(t,ci, label = 'Cl(intra)')
	plt.plot(t,(160- bl), label = 'Cl(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.xlabel('time (min)')
	plt.ylabel('Chloride Conc. (mM)')
	plt.legend(loc = 'right')

	plt.savefig(filename) # store local copy for later use
	plt.show()
	return filename

def vol_rat_calc(input_dict, t_on, t_off, t_end, volumes):
	# Wrapper function to pass dictionary parameters through solve_ivp
	def wrapper_fxn(t, y):
		return duct_model_system(t, y, cond)
	# Copy Input Dictionary so Initial Dictionary is Unaltered
	cond = copy.deepcopy(input_dict)

	model_results = []

	def reset_values():
		cond['bi'] = init_cond['bi']
		cond['bl'] = init_cond['bl']
		cond['ci'] = init_cond['ci']
		cond['ni'] = init_cond['ni']
		cond['gcftr'] = init_cond['gcftr']

	for volume in volumes:
		cond['vr'] = volume
		print(cond['vr'])
		# Period Before gcftr opens
		cond['gcftr'] = init_cond['gcftrbase']
		y0_0 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state0 = solve_ivp(wrapper_fxn, [0,t_on], y0_0)
		# Period where gcftr is open
		# Turn gcftr 'on' and re-initialize boundary conditions
		cond['gcftr'] = cond['gcftron']
		cond['bi'] = state0['y'][0][-1]
		cond['bl'] = state0['y'][1][-1]
		cond['ci'] = state0['y'][2][-1]
		cond['ni'] = state0['y'][3][-1]
		y0_1 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state1 = solve_ivp(wrapper_fxn, [t_on,t_off], y0_1)
		# Period where gcftr closes to end
		# Turn gcftr 'on' and re-initialize boundary conditions
		cond['gcftr'] = cond['gcftrbase']
		cond['bi'] = state1['y'][0][-1]
		cond['bl'] = state1['y'][1][-1]
		cond['ci'] = state1['y'][2][-1]
		cond['ni'] = state1['y'][3][-1]
		y0_2 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state2 = solve_ivp(wrapper_fxn, [t_off,t_end], y0_2)

		t_graph = np.concatenate([state0['t'], state1['t'], state2['t']])
		bi_graph = np.concatenate([state0['y'][0], state1['y'][0], state2['y'][0]])
		bl_graph = np.concatenate([state0['y'][1], state1['y'][1], state2['y'][1]])
		ci_graph = np.concatenate([state0['y'][2], state1['y'][2], state2['y'][2]])
		ni_graph = np.concatenate([state0['y'][3], state1['y'][3], state2['y'][3]])
		model_results.append([[t_graph, bi_graph, bl_graph, ci_graph, ni_graph], t_on, t_off, t_end, volume])
		reset_values()

	return model_results

def graph_volume_ratios(model_results, filename):
	for volume_option in model_results:
		t = volume_option[0][0]
		bl = volume_option[0][2]
		vol_label = 'VR ' + str(volume_option[-1])
		plt.plot(t, bl, label = vol_label)
	plt.xlabel('time (min)')
	plt.ylabel('HCO3- Conc. (mM)')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.title('Effect of different duct cell/lumen volume ratio on SS HCO3-')
	plt.savefig(filename, transparent=True) # store local copy for later use
	plt.show()
	return filename

def antiporters_calc(input_dict, t_on, t_off, t_end):
	# Wrapper function to pass dictionary parameters through solve_ivp
	def wrapper_fxn(t, y):
		return duct_model_system(t, y, cond)
	# Copy Input Dictionary so Initial Dictionary is Unaltered
	cond = copy.deepcopy(input_dict)
	model_results = []
	def reset_values():
		cond['bi'] = init_cond['bi']
		cond['bl'] = init_cond['bl']
		cond['ci'] = init_cond['ci']
		cond['ni'] = init_cond['ni']
		cond['gcftr'] = init_cond['gcftr']

	antiporter_options = [(True, True), (False, False), (True, False)]

	for option in antiporter_options:
		cond['ap_status'] = option[0]
		cond['apb_status'] = option[1]
		print(cond['ap_status'], cond['apb_status'])
		# Period Before gcftr opens
		cond['gcftr'] = init_cond['gcftrbase']
		y0_0 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state0 = solve_ivp(wrapper_fxn, [0,t_on], y0_0)
		# Period where gcftr is open
		# Turn gcftr 'on' and re-initialize boundary conditions
		cond['gcftr'] = cond['gcftron']
		cond['bi'] = state0['y'][0][-1]
		cond['bl'] = state0['y'][1][-1]
		cond['ci'] = state0['y'][2][-1]
		cond['ni'] = state0['y'][3][-1]
		y0_1 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state1 = solve_ivp(wrapper_fxn, [t_on,t_off], y0_1)
		# Period where gcftr closes to end
		# Turn gcftr 'on' and re-initialize boundary conditions
		cond['gcftr'] = cond['gcftrbase']
		cond['bi'] = state1['y'][0][-1]
		cond['bl'] = state1['y'][1][-1]
		cond['ci'] = state1['y'][2][-1]
		cond['ni'] = state1['y'][3][-1]
		y0_2 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
		state2 = solve_ivp(wrapper_fxn, [t_off,t_end], y0_2)

		t_graph = np.concatenate([state0['t'], state1['t'], state2['t']])
		bi_graph = np.concatenate([state0['y'][0], state1['y'][0], state2['y'][0]])
		bl_graph = np.concatenate([state0['y'][1], state1['y'][1], state2['y'][1]])
		ci_graph = np.concatenate([state0['y'][2], state1['y'][2], state2['y'][2]])
		ni_graph = np.concatenate([state0['y'][3], state1['y'][3], state2['y'][3]])
		model_results.append([[t_graph, bi_graph, bl_graph, ci_graph, ni_graph], t_on, t_off, t_end, option])
		reset_values()

	return model_results

def graph_antiporters(model_results, filename):
	for i in range(len(model_results)):
		t = model_results[i][0][0]
		bi = model_results[i][0][1]
		bl = model_results[i][0][2]
		ci = model_results[i][0][3]
		ni = model_results[i][0][4]
		t_on = model_results[i][1]
		t_off = model_results[i][2]
		t_end = model_results[i][3]
		option = model_results[i][4]

		# Three Subplots
		plt.subplot(3, 1, i+1)
		plt.plot(t,bi, 'r-', label = 'HCO3-(intra)')
		plt.plot(t,bl, 'b-', label = 'HCO3-(luminal)')
		plt.axvline(x=t_on, color = 'pink', label = 't_on')
		plt.axvline(x=t_off, color = 'purple', label = 't_off')
		plt.axvspan(t_on, t_off, alpha=0.1, color='red')
		plt.legend(loc = 'right')
		plt.ylim((0,155))
		plt.ylabel('Bicarb Conc. (mM)')
		if i+1 == 1:
			plt.title('Antiporters ON/OFF')

	plt.savefig(filename, transparent = True)
	plt.show()

	return filename



