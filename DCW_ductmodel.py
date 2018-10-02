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
vr = 0.1
bi0 = 15
buf = 0.1
chi = 1
bi = 15
ci = 60
bl = 32
cl = 160 - bl
gcftron = 1
gcftrbase = 0.00007
ek = -0.85
gk = 1
cap = 1
gnak = 3.125
np0 = 25
epump = -0.2
ionstr = 160
gnaleak = 0.4
jac = 0.25
rat = 0.25

def ductmodelsystem(state, t):
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
    jnak = gnak*(v-epump)*(ni/np0)
    jnaleak = gnaleak*(v-ena)
    flow = jlum*ionstr

    dbi = zeta*chi*(jbcftr+japl+japbl+buf*(bi0-bi)+2*jnbc)
    dbl = (jbl-jlum*bl)*zeta
    dci = jci*zeta
    dni = zeta*(jnbc-jnak-jnaleak)
    dgcftr = 0
    
    return [dbi, dbl, dci, dni, dgcftr]
    
t = np.linspace(0, 4000, 1000)
init_state = [bi, bl, ci, ni, gcftrbase] # 15, 32, 60, 28, 1 
state = odeint(ductmodelsystem, init_state, t)

plt.subplot(2, 1, 1)
plt.plot(t,state[:,0], 'r-', label = 'b_intra')
plt.plot(t,state[:,1], 'g-', label = 'b_luminal')
plt.legend()
plt.ylim((0,150))
plt.title('Duct Modeling Dif. Eq.')
plt.ylabel('Bicarb Conc. (mM)')

plt.subplot(2, 1, 2)
plt.plot(t,state[:,2], label = 'c_intra')
plt.plot(t,(160- state[:,1]), label = 'c_luminal')
plt.xlabel('time (min)')
plt.ylabel('Chloride Conc. (mM)')
plt.legend()
plt.ylim((0,150))
plt.show()


