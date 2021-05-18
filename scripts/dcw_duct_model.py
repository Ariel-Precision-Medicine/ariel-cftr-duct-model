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
			  'rat': 0.25, 'variant_adj': None, 'vr': 0.1, 'apb_status': False,
			  'ap_status': False, 'gcftr': 0.00007, 'smoke_adj': None,
			  'alcohol_adj': None}


def antiporter(ao,ai,bo,bi,ka,kb):
	'''
	Antiporter Function (ap(ao,ai,bo,bi,ka,kb) in original XPPAUT)

	Based on derived equation in Sohma 1996 ("A Mathematical Model of Pancreatic Ductal Epithelium").

	Assumptions for Binding Kinetics
	--------------------------------
	* antiporters had a single site which bound the transported ions
	* transported ions competed for the single binding site
	* that the antiporters did not cross the membrane in the absence of a bound ion
	* velocity constants for transport from the outside to inside and vice versa were the same
	* that the dissociation constants (ka,kb) for each ion at the intracellular and extracellular faces of the membrane were the same


	Parameters
	----------
	ao : float
		Molar concentration of ion A outside of cell membrane
	ai : float
		Molar concentration of ion A inside cell membrane
	bo : float
		Molar concentration of ion B outside of cell membrane
	bi : float
		Molar concentration of ion B inside of cell membrane
	ka : float
		Dissociation constant for ion A
	kb : float
		Dissociation constant for ion B

	Returns
	-------
	ratio : float
		Fraction of ionic flux moving through antiporter. Multiply with conductance to yield flux through antiporter.
	'''
    numerator = ao*bi-bo*ai
    denominator = (ka*kb*((1+ai/ka+bi/kb)*(ao/ka+bo/kb)+\
                          (1+ao/ka+bo/kb)*(ai/ka+bi/kb)))
    ratio = numerator / denominator
    return ratio


def eff_perm(xi,xo):
	'''
	Effective Permeability Function (g(xi,xo) in original).
	
	Parameters
	----------
	xi : float
		Molar concentration of ion X inside of the cell membrane
	xo : float
		Molar concentration of ion X outside of the cell membrane

	Returns
	-------
	float
		Coefficient found by linearizing the constant field equation around the equilibrium potential

	'''
    return (xi*xo*np.log(xi/xo)/(xi-xo))

def nernst_potential(a,b):
	'''
	Nernst Potential Function

	Parameters
	----------
	a : float
		Concentration inside cell
	b : float
		Concentration outside cell

	Returns
	-------
	float
		Nernst potential in volts (V)
	'''
    # Physical Constants
    ideal_gas = 8.31451
    faraday_cst = 96485
    body_temp = 310 #K
    return (ideal_gas*body_temp/faraday_cst)*np.log(a/b)

def duct_model_system(t, y, cond):
	'''
	System of ordinary differential equations to model ion flux of pancreatic ductal epithelial cells
	
	Parameters
	----------
	t : np.array
		Array of time series to solve differential equations along
	y : np.array
		Initial conditions of the system:
		* y[0] = intracellular bicarbonate concentration (HCO3-)
		* y[1] = luminal bicarbonate concentration (HCO3-)
		* y[2] = intracellular chloride concentration (Cl-)
		* y[3] = intracellular sodium concentration (Na+)
		* y[4] = state of CFTR channel (OPEN = 1, closed = ~0)
		Luminal chloride is linked to y[1] so it is not explicitly described.
		Similarly, intracellular potassium is linked to y[3] so it is not explicitly described.
	cond: dict
		Dictionary containing experimentally-derived parameters described in original Ermentrout paper.
		Open for updates over time as more data become available.

	Returns
	-------
	array
		Nested arrays to describe the change in ion concentration at each time step in the simulation.

	'''

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
	variant_adj = cond['variant_adj']
	smoke_adj = cond['smoke_adj']
	alcohol_adj = cond['alcohol_adj']

	if variant_adj == None:
		variant_adj = 1

	if smoke_adj == None:
		smoke_adj = 1

	if alcohol_adj == None:
		alcohol_adj = 1


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
	jbcftr = kbcf*(v-eb) * smoke_adj
	jccftr = (kccf*(v-ec) * variant_adj) * alcohol_adj

	# Antiporter Status (Open = True, Close = False)
	if ap_status:
		japl = antiporter(bl,bi,cl,ci,kbi,kcl)*gapl
	else:
		japl = 0
	if apb_status:
		japbl = antiporter(bb,bi,cb,ci,kbi,kcl)*gapbl
	else:
		japbl = 0

	jbl = (-jbcftr-japl)/vr+jac*rat
	jci = jccftr-japl-japbl
	jcl = ((-jccftr+japl)/vr+jac) 
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




