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
import numpy.random as rnd

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
ek = -0.085 
gk = 1
cap = 1
gnak = 3.125
np0 = 25
epump = -0.2
ionstr = 160
gnaleak = 0.4
jac = 0.025 
rat = 0.25

# cond_adj = 1 # adjust for functionality of different variants of CFTR, Cutting Paper

def ductmodelsystem(state, t, vr, ap_status, apb_status, g_bi, g_cl, cond_adj):
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
    jccftr = kccf*(v-ec) * cond_adj 
    japl = antiporter(bl,bi,cl,ci,kbi,kcl)*gapl*ap_status
    japbl = antiporter(bb,bi,cb,ci,kbi,kcl)*gapbl*apb_status
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
	apb_status = 1
	ap_status = 1
	g_bi = 0.2
	g_cl = 1
	cond_adj = 1

	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on) 
	state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,g_bi,g_cl,cond_adj,))
	on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
	t2 = np.linspace(t_on, t_off)
	state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,g_bi,g_cl,cond_adj,))
	off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
	t3 = np.linspace(t_off, 200000)
	state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi,g_cl,cond_adj,))

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
	apb_status = 1
	ap_status = 1
	cond_adj = 1

	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on)
	t2 = np.linspace(t_on, t_off)
	t3 = np.linspace(t_off, 200000)
	t = np.append(t1, t2)
	t = np.append(t,t3)
	t /= 20000

	i = 0
	# Repeat ODE simulation for each instance of volume ratio
	for vr in volume_ratios:
		state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status, g_bi, g_cl, cond_adj,))
		on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
		state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status, g_bi, g_cl, cond_adj,))
		off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
		state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
		state = np.vstack((state1, state2))
		state = np.vstack((state, state3))
		states[i] = state
		i += 1

	j = 0
	# Plot each line on the same graph
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

def graph_antiporters(init_state, filename):
	vr = 0.1
	states = [None, None, None]
	cond_adj = 1

	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on)
	t2 = np.linspace(t_on, t_off)
	t3 = np.linspace(t_off, 200000)
	t = np.append(t1, t2)
	t = np.append(t,t3)

	# Both antiporters on
	ap_status, apb_status = 1, 1
	state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
	state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
	state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	state = np.vstack((state1, state2))
	state = np.vstack((state, state3))

	states[0] = state

	# Both antiporters off
	ap_status, apb_status = 0, 0
	state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
	state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
	state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	state = np.vstack((state1, state2))
	state = np.vstack((state, state3))

	states[1] = state

	# APb off
	ap_status, apb_status = 0, 1
	state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
	state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
	state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi, g_cl, cond_adj,))
	state = np.vstack((state1, state2))
	state = np.vstack((state, state3))

	states[2] = state

	# Graph 3 subplots
	plt.subplot(3, 1, 1)
	plt.plot(t,states[0][:,2], 'r-', label = 'HCO3-(intra)')
	plt.plot(t,states[0][:,1], 'b-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.title('Antiporters ON/OFF')
	plt.ylabel('Bicarb Conc. (mM)')

	plt.subplot(3, 1, 2)
	plt.plot(t,states[1][:,2], 'r-', label = 'HCO3-(intra)')
	plt.plot(t,states[1][:,1], 'b-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.ylabel('Bicarb Conc. (mM)')

	plt.subplot(3, 1, 3)
	plt.plot(t,states[1][:,2], 'r-', label = 'HCO3-(intra)')
	plt.plot(t,states[1][:,1], 'b-', label = 'HCO3-(luminal)')
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.legend(loc = 'right')
	plt.ylim((0,155))
	plt.ylabel('Bicarb Conc. (mM)')

	plt.savefig(filename, transparent = True)
	plt.show()

	return

def graph_line_plots(init_state, filename):
	fig = plt.figure()

	vr = 0.1
	apb_status = 1
	ap_status = 1
	g_bi = [0.001, 0.2, 1]
	g_cl = [0.001, 0.01, 1]
	cond_adj = 1

	result_bl = []

	for bi_option in g_bi:
		for cl_option in g_cl:
			t_on, t_off = 20000, 120000
			t1 = np.linspace(0, t_on) 
			state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,bi_option,cl_option, cond_adj,))
			on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
			t2 = np.linspace(t_on, t_off)
			state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,bi_option,cl_option, cond_adj,))
			off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
			t3 = np.linspace(t_off, 200000)
			state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,bi_option,cl_option, cond_adj,))

			t = np.append(t1, t2)
			t = np.append(t,t3)
			state = np.vstack((state1, state2))
			state = np.vstack((state, state3))

			result_bl.append(state[-1,1])

	print(result_bl)

	# Process results collect from nested for loops
	# Place each 3rd item in appropriate group
	line_001 = result_bl[::3]
	line_02 = result_bl[1::3]
	line_1 = result_bl[2::3]

	print(line_001)
	print(line_02)
	print(line_1)




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

	plt.subplot(221)
	plt.plot(g_bi, line_001, label = str(g_bi[0]))
	plt.plot(g_bi, line_02, label = str(g_bi[1]))
	plt.plot(g_bi, line_1, label = str(g_bi[2]))
	plt.ylabel('Luminal Bicarb')
	plt.xlabel('gcftr (bi)')
	plt.legend(loc = 'right')


	plt.subplot(222)
	plt.imshow(rnd.random((100, 100)))


	plt.subplot(223)
	plt.imshow(rnd.random((100, 100)))


	plt.subplot(224)
	plt.imshow(rnd.random((100, 100)))

	plt.savefig(filename, transparent = True)
	plt.show()
	return

def visualize_conductance_effect(init_state, filename):
	vr = 0.1
	apb_status = 1
	ap_status = 1
	g_bi = 0.2
	g_cl = 1
	cond_adj = [0.054, 0.076, 0.151, 0.224, 0.563, 1]
	var_names = ['Q98R','F311L','F1099L', 'P5L', 'R31L', 'WT']
	states = [None] * len(cond_adj)

	t_on, t_off = 20000, 120000
	t1 = np.linspace(0, t_on)
	t2 = np.linspace(t_on, t_off)
	t3 = np.linspace(t_off, 200000)
	t = np.append(t1, t2)
	t = np.append(t,t3)
	t /= 2000
	t_on /= 2000
	t_off /= 2000

	i = 0
	# Repeat ODE for different conductance values
	for condition in cond_adj:
		state1 = odeint(ductmodelsystem, init_state, t1, args=(vr,ap_status,apb_status,g_bi,g_cl,condition,))
		on_state = [state1[-1,0],state1[-1,1],state1[-1,2],state1[-1,3], gcftron]
		state2 = odeint(ductmodelsystem, on_state, t2, args=(vr,ap_status,apb_status,g_bi,g_cl,condition,))
		off_state = [state2[-1,0],state2[-1,1],state2[-1,2],state2[-1,3], gcftrbase]
		state3 = odeint(ductmodelsystem, off_state, t3, args=(vr,ap_status,apb_status,g_bi,g_cl,condition,))
		state = np.vstack((state1, state2))
		state = np.vstack((state, state3))
		states[i] = state
		i += 1

	j = 0
	# Plot lines
	for state in states:
		plt.subplot(2, 1, 1)
		plt.plot(t,state[:,0], color = (1*cond_adj[j],0,0,1*cond_adj[j]), label = '(intra) ' + var_names[j] + '(' +str(round(cond_adj[j]*100)) +'%)' )
		plt.plot(t,state[:,1], color = (0,1*cond_adj[j],0,1*cond_adj[j]), label = '(luminal) ' +  var_names[j] + '(' +str(round(cond_adj[j]*100)) +'%)')

		plt.subplot(2, 1, 2)
		plt.plot(t,state[:,2], color = (0,0,1*cond_adj[j],1*cond_adj[j]), label = '(intra) '+  var_names[j] + '(' +str(round(cond_adj[j]*100)) +'%)')
		plt.plot(t,(160- state[:,1]), color = (1*cond_adj[j],1*cond_adj[j], 0,1*cond_adj[j]), label = '(luminal) '+  var_names[j] + '(' +str(round(cond_adj[j]*100)) +'%)')

		plt.ylim((0,155))
		j += 1

	plt.subplot(2, 1, 1)
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.axvspan(t[0], t_on, alpha=0.1, color='blue')
	plt.axvspan(t_off, t[-1], alpha=0.1, color='blue')
	plt.legend(loc = 'right', bbox_to_anchor= (1.1,0.5))
	plt.ylim((0,155))
	plt.title('Variant Conductivity Differences in Duct Modeling Dif. Eq. \n GCFTR ON in RED // GCFTR OFF in BLUE')
	plt.ylabel('Bicarb Conc. (mM)')

	plt.subplot(2, 1, 2)
	plt.axvline(x=t_on, color = 'pink', label = 't_on')
	plt.axvline(x=t_off, color = 'purple', label = 't_off')
	plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	plt.axvspan(t[0], t_on, alpha=0.1, color='blue')
	plt.axvspan(t_off, t[-1], alpha=0.1, color='blue')
	plt.ylabel('Chloride Conc. (mM)')
	plt.legend(loc = 'right', bbox_to_anchor=(1.1,0.5))
	plt.xlabel('time (min)')



	plt.savefig(filename, transparent=True) # store local copy for later use
	plt.show()
	return


### Not Working (debugging needed or unfinished) ###
#graph_antiporters([bi, bl, ci, ni, gcftrbase], 'antiporters.png')
#graph_line_plots([bi, bl, ci, ni, gcftrbase], 'quadrant.png')
#visualize_conductance_effect([bi, bl, ci, ni, gcftrbase], 'conductance_options.png')

### Working ###
#vol_rat_plot([bi, bl, ci, ni, gcftrbase], 'volumes.png')
#graph_HCO3_and_Cl([bi, bl, ci, ni, gcftrbase], 'smoke.png') # 15, 32, 60, 28, 1

'''
	smoking_penalty_rate = 0.15
	smoking_state = [bi*smoking_penalty_rate, bl*smoking_penalty_rate, ci, ni, gcftrbase]
	smoking_state1 = odeint(ductmodelsystem, smoking_state, t1)'''
	
