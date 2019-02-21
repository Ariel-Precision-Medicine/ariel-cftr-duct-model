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
	cond_adj = cond['cond_adj']

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

