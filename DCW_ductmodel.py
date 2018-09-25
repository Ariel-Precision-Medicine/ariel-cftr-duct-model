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
    return (ideal_gas*body_temp/faraday_cst)*log(a,b)

# Voltage-Dependent Flux Fxn
def voltage_flux(perm, voltage_potential, nernst):
    return (perm*(voltage_potential-nernst))

# Conductances
g_nbc = 2
g_apl = 0.25
g_apbl = 0.005

# bicarb9.ode
g_bi = 0.2
g_cl = 1

# Time scale adjustment factor (multiplication of all permeabilities)
zeta = 0.05 # (lowercase xi in Eq. 21,22,23,24)

# Antiporter Constants
k_bi = 1
k_cl = 10

# Basolateral Concentrations
nb = 140 # sodium basolateral
bb = 22 # bicarbonate basolateral
cb = 130 # chloride basolateral

# Intracellular sodium (used to turn Na/K pump on/off)
ni = 14 # units?

# Volume Ratio (10:1), rho is ClBi permeability
volume_ratio = 0.1
buf = 0.1
chi = 1

# Initial Intracellular concentrations
intra_cl = 60 # intracellular chloride
luminal_bi = 32 # luminal bicarb
b_i_0 = 15 # intracellular bicarb
na_i_0 = 25 # intracellular sodium

# Osmotic Balance (160mM)
luminal_cl = 160 - luminal_bi
g_CFTR_on = 1
g_CFTR_base = 7*10**(-5)
ek = -0.085
gk = 1

# Membrane Capacitance
mem_capacitance = 1

# Voltage via Nernst Potential
eb = nernst_potential(b_i_0, luminal_bi)
enbc = nernst_potential((b_i_0**2*ni), (bb**2*nb))
ec = nernst_potential(intra_cl, luminal_cl)
ena = nernst_potential(nb, ni)

gnaleak = 0.4

# Acinar Variables
jac = 0.025 # Jacobian??? Flux??
acinar_ratio = 0.25 # Maintain ratio of 4:1 Cl:Bi (original "rat")

# Ionic Strength
ionic_strength = 160
g_nak = 3.125
e_pump = -0.2

# effective CFTR permeability
knbc = g_nbc #  What is this for?


###############################################################################
# Differential Equations for the Concentrations
###############################################################################

# Times from Figure 3 graphs
time_on = 1 # gate opens (seemingly arbitrary)
time_off = 6 # gate closes (seemingly arbitrary)
time_start = 0 
time_end = 10

# Initial Intracellular concentrations
### Add basolateral concentrations to act as initial conditions
intra_cl = 60 # intracellular chloride
luminal_bi = 32 # luminal bicarb
b_i_0 = 15 # intracellular bicarb
na_i_0 = 25 # intracellular sodium
gcftr_0 = 0 # changed to 1 to "turn on" according to paper


def flux(gcftr, bi, t):
    # Effective CFTR Permeability
    k_bi_CFTR = eff_perm(bi, luminal_bi)*gcftr*g_bi
    k_cl_CFTR = eff_perm(intra_cl, luminal_cl)*gcftr*g_cl
    # Voltage Potential
    voltage_potential = (knbc*enbc + k_bi_CFTR*eb + k_cl_CFTR*ec + gk*ek + \
                     gnaleak*ena) / (knbc + k_bi_CFTR + k_cl_CFTR + gk)
    # All fluxes
    J_nbc = voltage_flux(knbc, voltage_potential, enbc)
    J_ccftr = voltage_flux(k_cl_CFTR, voltage_potential, ec)
    J_bcftr = voltage_flux(k_bi_CFTR, voltage_potential, eb)
    J_apl = antiporter(luminal_bi,bi,luminal_cl,intra_cl,k_bi,k_cl) * g_apl
    J_apbl = antiporter(bb,bi,cb,intra_cl,k_bi,k_cl) * g_apbl
    J_intra_cl = J_ccftr - J_apl - J_apbl
    J_luminal_bi = (-J_bcftr-J_apl)/volume_ratio + jac * acinar_ratio
    J_intra_cl = J_ccftr - J_apl - J_apbl
    J_luminal_cl = (-J_ccftr + J_apl)/volume_ratio + jac
    J_luminal = (J_luminal_cl + J_luminal_bi) / ionic_strength
    J_nak = g_nak * (voltage_potential - ena)
    J_na_leak = gnaleak * (voltage_potential - ena)
    
    # Differential Equations to be integrated & graphed
    
    # Luminal Cl remains at 160
    # form: dxdt = 0 (straight line, no integration needed)
    dcl_lumdt = zeta * J_luminal_cl  #mM (assumed to be constant) 
    ### Sum of luminal chloride and bicarbonate
    
    # Intracellular Cl
    # form: dxdt = a*b
    # dCl(intra)/dt = zeta * J_intra_cl
    dcl_intradt = zeta * J_intra_cl
    
    # Luminal Bicarb
    # form: dxdt = a*(b-c*x)
    # dBi(luminal)/dt = zeta * (J_luminal_bi - J_luminal * luminal_bi)
    dbi_lumdt = zeta * (J_luminal_bi - J_luminal * luminal_bi)
    
    # Intracellular Bicarb
    # form: dxdt = a*b*(c + d + e + f*(x_i-x) + g)
    # dBi(intra)/dt = zeta * chi *
    #                       (J_bcftr + J_apl + J_apbl + buf * (bi_0-bi) + 2*J_nbc)
    dbi_intradt = zeta * chi * (J_bcftr + J_apl + J_apbl + buf * (b_i_0-bi) + 2*J_nbc)
    
    # Intracellular Sodium
    # form: dxdt = a*(b - c - d)
    # dNa(intra)/dt = zeta * (J_nbc - J_nak - J_na_leak)
    dna_intradt = zeta * (J_nbc - J_nak - J_na_leak)
    
    return dcl_lumdt,dcl_intradt,dbi_lumdt,dbi_intradt,dna_intradt

## time points
#t = np.linspace(0,600)
#
## solve ODE
#Y0 = [0,0,0,0,0]
#y = odeint(flux, t, Y0)
#
## plot results
#plt.plot(t,y[0])
#plt.xlabel('time')
#plt.ylabel('y(t)')
#plt.show()
#
    
# Basolateral Concentrations
nb = 140 # sodium basolateral
bb = 22 # bicarbonate basolateral
cb = 130 # chloride basolateral
intra_cl = 60 # intracellular chloride
luminal_bi = 32 # luminal bicarb
b_i_0 = 15 # intracellular bicarb
na_i_0 = 25 # intracellular sodium
gcftr_0 = 0 # changed to 1 to "turn on" according to paper

t = np.linspace(0,600)
nb_graph = np.linspace(nb,nb)
bb_graph = np.linspace(bb,bb)
cb_graph = np.linspace(cb,cb)

plt.plot(t,nb_graph,label='basolateral Na')
plt.plot(t,bb_graph,label='basolateral Bicarb')
plt.plot(t,cb_graph,label='basolateral Cl')
plt.xlabel('time')
plt.ylabel('Concentration (mM)')
plt.legend(loc='upper right')
plt.show()
