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
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

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

# Constants & Initial Conditions
g_bi = 0.2
g_cl = 1
zeta = 0.05
kbi = 1
kcl = 10
gnbc = 2
gapl = 0.25
gapbl = 0.005
nb = 140
bb = 22
cb = 130
ni = 14
bi0 = 15
buf = 0.1
chi = 1
bi = 15
ci = 60
bl = 32
cl = 160 - bl
gcftron = 1
gcftrbase = 0.00007
ek = -0.085 # should be -0.085, but system behaves better with this number
gk = 1
cap = 1
gnak = 3.125
np0 = 25
epump = -0.2
ionstr = 160
gnaleak = 0.4
jac = 0.025 # should be 0.025
rat = 0.25

def ductmodelsystem(state, t, vr):
    bi, bl, ci, ni, gcftr = state
    cl = 160 - bl
    eb = nernst_potential(bi, bl)
    enbc = nernst_potential((bi**2*ni), (bb**2*nb))
    ec = nernst_potential(ci, cl)
    ena = nernst_potential(nb, ni)
    kccf = eff_perm(ci,cl)*gcftr*g_cl
    kbcf = eff_perm(bi, bl)*gcftr*g_bi
    knbc = gnbc
    v = (knbc*enbc+kbcf*eb+kccf*ec+gk*ek+gnaleak*ena)/(knbc+kbcf+kccf+gk)
    jnbc = knbc*(v-enbc)
    jbcftr = kbcf*(v-eb)
    jccftr = kccf*(v-ec)
    japl = antiporter(bl,bi,cl,ci,kbi,kcl)*gapl
    japbl = antiporter(bb,bi,cb,ci,kbi,kcl)*gapbl
    jbl = (-jbcftr-japl)/vr+jac*rat
    jci = jccftr-japl-japbl
    jcl = (-jccftr+japl)/vr+jac
    jlum = (jcl+jbl)/ionstr
    jnak = gnak*(v-epump)*(ni/np0)**3
    jnaleak = gnaleak*(v-ena)
    flow = jlum*ionstr

    dbi = zeta*chi*(jbcftr+japl+japbl+buf*(bi0-bi)+2*jnbc)
    dbl = (jbl-jlum*bl)*zeta
    dci = jci*zeta
    dni = zeta*(jnbc-jnak-jnaleak)
    dgcftr = 0
    return [dbi, dbl, dci, dni, dgcftr]

def graph_HCO3_and_Cl(init_state, filename):
	vr = 0.1
	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on) 
	state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,))
	on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
	t2 = np.linspace(t_on, t_off)
	state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,))
	off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
	t3 = np.linspace(t_off, 200000)
	state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,))

	t = np.append(t1, t2)
	t = np.append(t,t3)
	state = np.vstack((state1, state2))
	state = np.vstack((state, state3))

	plt.subplot(2, 1, 1)
	plt.plot(t,state[:,0], 'r-', label = 'HCO3-(intra)')
	plt.plot(t,state[:,1], 'g-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.title('Duct Modeling Dif. Eq. \n GCFTR ON in RED')
	plt.ylabel('Bicarb Conc. (mM)')

	plt.subplot(2, 1, 2)
	plt.plot(t,state[:,2], label = 'Cl(intra)')
	plt.plot(t,(160- state[:,1]), label = 'Cl(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.xlabel('time (min)')
	plt.ylabel('Chloride Conc. (mM)')
	plt.legend(loc = 'right')
	plt.ylim((0,155))

	plt.savefig(filename, transparent=True) # store local copy for later use
	plt.show()
	return

def vol_rat_plot(init_state, filename):
	volume_ratios = [10, 1, 0.1, 0.04]
	states = [None] * len(volume_ratios)

	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on)
	t2 = np.linspace(t_on, t_off)
	t3 = np.linspace(t_off, 200000)
	t = np.append(t1, t2)
	t = np.append(t,t3)
	t /= 20000

	i = 0
	for vr in volume_ratios:
		state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,))
		on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
		state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,))
		off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
		state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,))
		state = np.vstack((state1, state2))
		state = np.vstack((state, state3))
		states[i] = state
		i += 1

	j = 0
	for state in states:
		plt.plot(t, state[:, 1], label='VR'+str(volume_ratios[j]))
		j += 1

	plt.xlabel('time (min)')
	plt.ylabel('HCO3- Conc. (mM)')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.title('Effect of different duct cell/lumen VR on SS HCO3-')
	plt.savefig(filename, transparent=True) # store local copy for later use
	plt.show()
	return

vol_rat_plot([bi, bl, ci, ni, gcftrbase], 'volumes.png')
graph_HCO3_and_Cl([bi, bl, ci, ni, gcftrbase], 'nosmoke.png') # 15, 32, 60, 28, 1

'''
	smoking_penalty_rate = 0.15
	smoking_state = [bi*smoking_penalty_rate, bl*smoking_penalty_rate, ci, ni, gcftrbase]
	smoking_state1 = odeint(ductmodelsystem, smoking_state, t1)'''
	
