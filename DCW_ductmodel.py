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

#######
#    # Constant Vars
#    bi = 15 
#    luminal_bi = 32
#    g_bi = 0.2
#    intra_cl = 60 
#    luminal_cl = 160 - luminal_bi
#    g_cl = 1
#    g_nbc = 2
#    knbc = g_nbc
#    ni = 14
#    nb = 140 # sodium basolateral
#    bb = 22 # bicarbonate basolateral 
#    ek = -0.085
#    gk = 1
#    gnaleak = 0.4
#    enbc = nernst_potential((bi**2*ni), (bb**2*nb))
#    eb = nernst_potential(bi, luminal_bi)
#    ec = nernst_potential(intra_cl, luminal_cl)
#    ena = nernst_potential(nb, ni)
#    volume_ratio = 0.1
#    g_apl = 0.25
#    k_bi = 1
#    k_cl = 10
#    jac = 0.025
#    zeta = 0.05
#####
# Bicarb9 ODE
g_bi = 0.2
g_cl = 1

# Time scale adjustment factor (multiplication of all permeabilities)
zeta = 0.05

# Antiporter Fxn (ap(ao,aibo,bi,ka,kb) in original)
def antiporter(ao,ai,bo,bi,ka,kb):
    numerator = ao*bi-bo*ai
    denominator = (ka*kb*((1+ai/ka+bi/kb)*(ao/ka+bo/kb)+\
                          (1+ao/ka+bo/kb)*(ai/ka+bi/kb)))
    return (numerator/denominator)

# Antiporter Constants
k_bi = 1
k_cl = 10

# Effective Permeability Fxn (g(xi,xo) in original)
# Linearization of the Constant Field Eqn
def eff_perm(xi,xo):
    return (xi*xo*log(xi/xo)/(xi-xo))

# Conductances
g_nbc = 2
g_apl = 0.25
g_apbl = 0.005

# Basolateral Concentrations
nb = 140 # sodium basolateral
bb = 22 # bicarbonate basolateral
cb = 130 # chloride basolateral

# Intracellular sodium (used to turn Na/K pump on/off)
ni = 14 # units?

# Volume Ratio (10:1), rho is ClBi permeability
volume_ratio = 0.1
bi_0 = 15
buf = 0.1
chi = 1


# Intracellular concentrations
intra_cl = 60 # original "ci"
luminal_bi = 32 # original "bl"

# Osmotic Balance (160mM)
luminal_cl = 160 - luminal_bi
g_CFTR_on = 1
g_CFTR_base = 7*10**(-5)
ek = -0.085
gk = 1

# Membrane Capacitance
mem_capacitance = 1


# Nernst Potential Fxn
def nernst_potential(a,b):
    # Physical Constants
    ideal_gas = 8.31451
    faraday_cst = 96485
    body_temp = 310 #K
    return (ideal_gas*body_temp/faraday_cst)*log(a,b)

eb = nernst_potential(bi_0, luminal_bi)
enbc = nernst_potential((bi_0**2*ni), (bb**2*nb))
ec = nernst_potential(intra_cl, luminal_cl)
ena = nernst_potential(nb, ni)

# Effective CFTR Permeability
#k_cl_CFTR = eff_perm(intra_cl, luminal_cl)*gcftr*g_cl
#k_bi_CFTR = eff_perm(bi, luminal_bi)*gcftr*g_bi
knbc = g_nbc #  What is this for?

gnaleak = 0.4

# Voltage potential (slaved to concentrations)
#voltage_potential = (knbc*enbc + k_bi_CFTR*eb + k_cl_CFTR*ec + gk*ek + \
#                     gnaleak*ena) / (knbc + k_bi_CFTR + k_cl_CFTR + gk)

# Voltage-Dependent Flux Fxn
def voltage_flux(perm, voltage_potential, nernst):
    return (perm*(voltage_potential-nernst))

#J_nbc = voltage_flux(knbc, voltage_potential, enbc)
#J_bcftr = voltage_flux(k_bi_CFTR, voltage_potential, eb)
#J_ccftr = voltage_flux(k_cl_CFTR, voltage_potential, ec)

# Antiporter Flux
J_apl = antiporter(luminal_bi,bi_0,luminal_cl,intra_cl,k_bi,k_cl) * g_apl
J_apbl = antiporter(bb,bi_0,cb,intra_cl,k_bi,k_cl) * g_apbl

# Acinar Variables
jac = 0.025 # Jacobian???
acinar_ratio = 0.25 # Maintain ratio of 4:1 Cl:Bi (original "rat")

# Total Ionic Fluxes
#J_luminal_bi = (-J_bcftr-J_apl)/volume_ratio + jac * acinar_ratio
#J_intra_cl = J_ccftr - J_apl - J_apbl
#J_luminal_cl = (-J_ccftr + J_apl)/volume_ratio + jac

# Ionic Strength
ionic_strength = 160

# Total Scaled Luminal Flux
#J_luminal = (J_luminal_cl + J_luminal_bi) / ionic_strength

g_nak = 3.125
e_pump = -0.2
np_0 = 25

#J_nak = g_nak * (voltage_potential - ena)
#J_na_leak = gnaleak * (voltage_potential - ena)
#
## Flow
#flow = J_luminal * ionic_strength

# Differential Equations for the Concentrations

# Use Scipy ODE here
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html

time_on = 1 # gate opens (seemingly arbitrary)
time_off = 6 # gate closes (seemingly arbitrary)

time_start = 0 
time_end = 10




def flux(gcftr, bi, t):
    
    k_bi_CFTR = eff_perm(bi, luminal_bi)*gcftr*g_bi
    k_cl_CFTR = eff_perm(intra_cl, luminal_cl)*gcftr*g_cl
    
    voltage_potential = (knbc*enbc + k_bi_CFTR*eb + k_cl_CFTR*ec + gk*ek + \
                     gnaleak*ena) / (knbc + k_bi_CFTR + k_cl_CFTR + gk)
    
    J_nbc = voltage_flux(knbc, voltage_potential, enbc)
    J_ccftr = voltage_flux(k_cl_CFTR, voltage_potential, ec)
    J_bcftr = voltage_flux(k_bi_CFTR, voltage_potential, eb)
    J_apl = antiporter(luminal_bi,bi,luminal_cl,intra_cl,k_bi,k_cl) * g_apl
    J_intra_cl = J_ccftr - J_apl - J_apbl
    J_luminal_bi = (-J_bcftr-J_apl)/volume_ratio + jac * acinar_ratio
    J_intra_cl = J_ccftr - J_apl - J_apbl
    J_luminal_cl = (-J_ccftr + J_apl)/volume_ratio + jac
    J_luminal = (J_luminal_cl + J_luminal_bi) / ionic_strength
    J_nak = g_nak * (voltage_potential - ena)
    J_na_leak = gnaleak * (voltage_potential - ena)
    ###
    
    # Luminal Cl remains at 160
    # form: dxdt = 0
    cl_lum = 160 #mM (assumed to be constant)
    
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
    dbi_intradt = zeta * chi * (J_bcftr + J_apl + J_apbl + buf * (bi_0-bi) + 2*J_nbc)
    
    # Intracellular Sodium
    # form: dxdt = a*(b - c - d)
    # dNa(intra)/dt = zeta * (J_nbc - J_nak - J_na_leak)
    dna_intradt = zeta * (J_nbc - J_nak - J_na_leak)
    
    return cl_lum,dcl_intradt,dbi_lumdt,dbi_intradt,dna_intradt

# time points
t = np.linspace(0,600)

# solve ODE
y = odeint(flux, 0, t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()


