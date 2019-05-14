# duct_model_bokeh_plotting.py

from math import log
from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rnd
import copy
from DCW_duct_model import duct_model_system, init_cond

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column, row


def cftr_calc_HCO3_Cl(input_dict, t_on, t_off, t_end):
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

def graph_CFTR_open_close(model_results, filename, title):
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
				  x_axis_label = 'time (min',
				  y_axis_label = 'Bicarb Conc. (mM)')

	plot_bicarb.vbar(x=np.average([t_on, t_off]), bottom=0, top=np.max(bl), 
         color='gray', width=t_off-t_on, 
         legend='CFTR Channel Open', fill_alpha = 0.1)

	plot_bicarb.line(t, bl, legend = 'Luminal Bicarbonate',
					 line_width = 3,
					 line_color = 'red')

	plot_bicarb.line(t, bi, legend = 'Intracellular Bicarbonate',
					 line_width = 3,
					 line_color = 'orange')

	# Create new plot with title and axis
	plot_chloride = figure(title = title,
				  x_axis_label = 'time (min',
				  y_axis_label = 'Chloride Conc. (mM)')

	plot_chloride.vbar(x=np.average([t_on, t_off]), bottom=0, top=np.max(bl), 
         color='gray', width=t_off-t_on, 
         legend='CFTR Channel Open', fill_alpha = 0.1)

	plot_chloride.line(t, cl, legend = 'Luminal Chloride',
					 line_width = 3,
					 line_color = 'navy')

	plot_chloride.line(t, ci, legend = 'Intracellular Chloride',
					 line_width = 3,
					 line_color = 'blue')

	

	# plot.line(t, cl, legend = 'L')
	# plot.line(t, ci)

	show(column(plot_bicarb, plot_chloride))

	# # Upper Plot
	# plt.subplot(2, 1, 1)
	# plt.plot(t,bi, 'r-', label = 'HCO3-(intra)')
	# plt.plot(t,bl, 'g-', label = 'HCO3-(luminal)')
	# plt.axvline(x=t_on, color = 'pink', label = 'CFTR open')
	# plt.axvline(x=t_off, color = 'purple', label = 'CFTR closed')
	# plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	# plt.legend(loc = 'right')
	# plt.title('Duct Modeling Dif. Eq. \n GCFTR OPEN in RED')
	# plt.ylabel('Bicarb Conc. (mM)')

	# # Lower Plot
	# plt.subplot(2, 1, 2)
	# plt.plot(t,ci, label = 'Cl(intra)')
	# plt.plot(t,(160-bl), label = 'Cl(luminal)')
	# plt.axvline(x=t_on, color = 'pink', label = 'CFTR open')
	# plt.axvline(x=t_off, color = 'purple', label = 'CFTR closed')
	# plt.axvspan(t_on, t_off, alpha=0.1, color='red')
	# plt.xlabel('c')
	# plt.ylabel('Chloride Conc. (mM)')
	# plt.legend(loc = 'right')

	# plt.savefig(filename) # store local copy for later use
	# plt.show()
	return filename

graph_CFTR_open_close(cftr_calc_HCO3_Cl(init_cond, 20000, 120000, 200000),'CFTR_plot', 'Duct Modeling Differential Equation Analysis')
