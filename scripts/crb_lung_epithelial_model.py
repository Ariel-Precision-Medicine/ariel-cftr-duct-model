# crb_lung_epithelial_model.py

"""
Created July 22nd 2019

Cameron Breze for Ariel Precision Medicine

Second cellular mechanistic model to highlight patterns of ion
transport in lung epithelial cells that contain CFTR. As with the
first model (dcw_duct_model.py), CFTR variants are incorporated to
show how a change in one functional piece influences the rest of the
model on the whole.

"""

from math import log
from scipy.integrate import odeint, solve_ivp
import numpy as np
import numpy.random as rnd
import copy
from dcw_duct_model import antiporter, eff_perm, nernst_potential


