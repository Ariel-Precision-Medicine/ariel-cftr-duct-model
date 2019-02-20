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

# Antiporter Fxn (ap(ao,aibo,bi,ka,kb) in original)
def antiporter(ao,ai,bo,bi,ka,kb):
    numerator = ao*bi-bo*ai
    denominator = (ka*kb*((1+ai/ka+bi/kb)*(ao/ka+bo/kb)+\
                          (1+ao/ka+bo/kb)*(ai/ka+bi/kb)))
    return (numerator/denominator)

# Effective Permeability Fxn (g(xi,xo) in original)
# Linearization of the Constant Field Eqn
def eff_perm(xi,xo):
    return (xi*xo*log(xi/xo)/(xi-xo))

# Nernst Potential Fxn
def nernst_potential(a,b):
    # Physical Constants
    ideal_gas = 8.31451
    faraday_cst = 96485
    body_temp = 310 #K
    return (ideal_gas*body_temp/faraday_cst)*log(a/b)

def duct_model_system(t, y, cond):
	# Unpack variables to be integrated
	bi = y[0]
	bl = y[1]
	ci = y[2]
	ni = y[3]
	gcftr = y[4]
	# Unpack constants
	cl = 160 - bl
	nb = cond['nb']
	bb = cond['bb']
	g_cl = cond['g_cl']
	g_bi = cond['g_bi']
	gnbc = cond['gnbc']
	zeta = cond['zeta']
	kbi = cond['kbi']
	kcl = cond['kcl']
	gnbc = cond['gnbc']
	gapl = cond['gapl']
	gapbl = cond['gapbl']
	cb = cond['cb']
	gcftron = cond['gcftron']
	gcftrbase = cond['gcftrbase']
	ek = cond['ek']
	gk = cond['gk']
	cap = cond['cap']
	gnak = cond['gnak']
	np0 = cond['np0']
	epump = cond['epump']
	ionstr = cond['ionstr']
	gnaleak = cond['gnaleak']
	jac = cond['jac']
	rat = cond['rat']
	vr = cond['vr']
	apb_status = cond['apb_status']
	ap_status = cond['ap_status']
	chi = cond['chi']
	buf = cond['buf']
	bi0 = cond['bi0']

	# Unpack Added Variables
	cond_adj = init_cond['cond_adj']

	# Nernst Potentials
	eb = nernst_potential(bi, bl)
	enbc = nernst_potential((bi**2*ni), (bb**2*nb))
	ec = nernst_potential(ci, cl)
	ena = nernst_potential(nb, ni)

	# Permeabilities
	kccf = eff_perm(ci,cl)*gcftr*g_cl 
	kbcf = eff_perm(bi, bl)*gcftr*g_bi
	knbc = gnbc

	# Voltage Potential
	v = (knbc*enbc+kbcf*eb+kccf*ec+gk*ek+gnaleak*ena)/(knbc+kbcf+kccf+gk)

	# Flux Calculations
	jnbc = knbc*(v-enbc)
	jbcftr = kbcf*(v-eb)
	jccftr = kccf*(v-ec) * cond_adj

	# Antiporter Status (Open = True, Close = False)
	if ap_status:
		japl = antiporter(bl,bi,cl,ci,kbi,kcl)*gapl
	else:
		japl = 0
	if apb_status:
		japbl = antiporter(bb,bi,cb,ci,kbi,kcl)*gapbl*apb_status
	else:
		japbl = 0

	jbl = (-jbcftr-japl)/vr+jac*rat
	jci = jccftr-japl-japbl
	jcl = (-jccftr+japl)/vr+jac
	jlum = (jcl+jbl)/ionstr
	jnak = gnak*(v-epump)*(ni/np0)**3
	jnaleak = gnaleak*(v-ena)
	flow = jlum*ionstr

	# Differential Equations
	dbi = zeta*chi*(jbcftr+japl+japbl+buf*(bi0-bi)+2*jnbc)
	dbl = (jbl-jlum*bl)*zeta
	dci = jci*zeta
	dni = zeta*(jnbc-jnak-jnaleak)
	dgcftr = 0
	return [dbi, dbl, dci, dni, dgcftr]


def cftr_calc_HCO3_Cl(input_dict, t_on, t_off, t_end, y0_0):
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

	plt.subplot(2, 1, 1)
	plt.plot(t,bi, 'r-', label = 'HCO3-(intra)')
	plt.plot(t,bl, 'g-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.title('Duct Modeling Dif. Eq. \n GCFTR OPEN in RED')
	plt.ylabel('Bicarb Conc. (mM)')

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

y0_0 = [init_cond['bi'], init_cond['bl'], init_cond['ci'], init_cond['ni'], init_cond['gcftr']]
graph_CFTR(cftr_calc_HCO3_Cl(init_cond, 20000, 120000, 200000, y0_0), 'testmodeling')













