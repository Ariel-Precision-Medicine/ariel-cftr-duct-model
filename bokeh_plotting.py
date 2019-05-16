# bokeh_plotting.py

from math import log
from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rnd
import copy
from DCW_duct_model import duct_model_system, init_cond

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column, row
from bokeh.models import Legend


def run_model_CFTR(input_dict, t_on, t_off, t_end):
	times = {'t_on': t_on, 't_off': t_off, 't_end': t_end}

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
	graphing_strings = ['bi', 'bl', 'ci', 'ni']
	for variable in graphing_strings:
		cond[variable] = state0['y'][graphing_strings.index(variable)][-1]
	
	y0_1 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
	state1 = solve_ivp(wrapper_fxn, [t_on,t_off], y0_1)

	# Period where gcftr closes to end
	# Turn gcftr 'on' and re-initialize boundary conditions
	cond['gcftr'] = cond['gcftrbase']
	for variable in graphing_strings:
		cond[variable] = state1['y'][graphing_strings.index(variable)][-1]

	y0_2 = [cond['bi'], cond['bl'], cond['ci'], cond['ni'], cond['gcftr']]
	state2 = solve_ivp(wrapper_fxn, [t_off,t_end], y0_2)

	graphing_dict = dict()
	graphing_dict['time'] = np.concatenate([state0['t'], state1['t'], state2['t']])
	for variable in graphing_strings:
		graphing_dict[variable] = np.concatenate([state0['y'][graphing_strings.index(variable)], 
									state1['y'][graphing_strings.index(variable)],
									state2['y'][graphing_strings.index(variable)]])
	return [graphing_dict, times, graphing_strings]

def graph_CFTR(model_results, filename, title):
	# Unpack variables
	graphing_dict, times, graphing_strings = model_results[0], model_results[1], model_results[2]
	time_adj = 20000
	t = graphing_dict['time'] / time_adj
	bi = graphing_dict['bi']
	bl = graphing_dict['bl']
	ci = graphing_dict['ci']
	ni = graphing_dict['ni']
	cl = 160 - bl # from DCW model assumptions
	t_on = times['t_on'] / time_adj
	t_off = times['t_off'] / time_adj
	t_end = times['t_end'] / time_adj

	# Create output file
	output_file(filename + '.html')
	# Create new plot with title and axis
	plot_bicarb = figure(title = title,
				  x_axis_label = 'time (min)',
				  y_axis_label = 'Bicarb Conc. (mM)',
				  y_range = (0, 150))

	channel_open = plot_bicarb.vbar(x=np.average([t_on, t_off]), bottom=0, top=np.max(bl), 
							         color='gray', width=t_off-t_on, fill_alpha = 0.1)

	lum_bicarb = plot_bicarb.line(t, bl,
										 line_width = 3,
										 line_color = '#34344A')

	intra_bicarb = plot_bicarb.line(t, bi,
											 line_width = 3,
											 line_color = '#7FE0CB')

	peak_bicarb = plot_bicarb.circle(t_off, np.max(bl),
														   size=20, 
														   color='red', 
														   alpha=0.5)

	secretion_text = 'Peak Secretion '+str(np.max(bl).round(0))+' mM'

	legend = Legend(items=[
			    	("Luminal Bicarbonate"   , [lum_bicarb]),
				    ("Intracellular Bicarbonate" , [intra_bicarb]),
				    (secretion_text , [peak_bicarb]),
				    ("CFTR Channel Open", [channel_open])
							], location="center")

	plot_bicarb.add_layout(legend, 'right')


	# Create new plot with title and axis
	plot_chloride = figure(title = title,
				  x_axis_label = 'time (min)',
				  y_axis_label = 'Chloride Conc. (mM)',
				  y_range = (0, 150))

	plot_chloride.vbar(x=np.average([t_on, t_off]), bottom=0, top=np.max(cl), 
         color='gray', width=t_off-t_on, 
         legend='CFTR Channel Open', fill_alpha = 0.1)

	plot_chloride.line(t, cl, legend = 'Luminal Chloride',
					 line_width = 3,
					 line_color = '#34344A')

	plot_chloride.line(t, ci, legend = 'Intracellular Chloride',
					 line_width = 3,
					 line_color = '#7FE0CB')

	plot_chloride.background_fill_color = '#F4F1E1'
	plot_chloride.background_fill_alpha = 0.5

	return plot_bicarb, plot_chloride

graph_CFTR(run_model_CFTR(init_cond, 20000, 120000, 200000),'CFTR_plot', 'Duct Modeling Differential Equation Analysis')