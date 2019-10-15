# beagle_text_processing.py

'''
Cameron Breze for Ariel Precision Medicine
Purpose: Processed Phased CFTR data from Beagle Output

Dependencies:
-pandas
-openpyxl
-myvariant

'''

import pandas as pd
import myvariant
import string
import copy
import math
import statistics
import os
from oop_duct_model import Duct_Cell
from bokeh.plotting import show, figure, save
from bokeh_plotting import run_model_CFTR
from dcw_duct_model import init_cond
from bokeh.models import Legend


def appendColumnHeaders(colTextFile,dataTextFile):
	# Read column text file and create df to hold information
	columns = pd.read_csv('inputs/' + colTextFile, sep='\t')
	# Pass columns alongside data set to merge into one df
	df = pd.read_csv('inputs/' + dataTextFile, sep='\t', names = columns)
	return df

def exportExcel(filename, dataframe):
	# Generate record of dataframe joining
	dataframe.to_excel(filename)
	return

def gatherPatientIDs(dataframe):
	# Gather columns from dataframe
	columns = list(dataframe.columns)
	# Index last location of common column headers
	start = columns.index('FORMAT') + 1
	# Generate list with all patient IDs
	patientIDs = columns[start:]
	return patientIDs

def processPhasedData(phasingPattern):
	# Input format '#|#', returns integer for L/R chromosome
	if phasingPattern == '0|0':
		# ignore null input, bug protection
		pass
	else:
		# separate the phased data by chromosome and change string to int
		splitPattern = phasingPattern.split("|")
		leftChromosome, rightChromosome = splitPattern[0], splitPattern[1]
	return leftChromosome, rightChromosome

def gatherPhasedVariants(patientID, df):
	# Returns dictionary of format below where each variant is added to
	# the appropriate phasing list
	phasedVariants = {'ID': patientID, 'ChrL': [], 'ChrR': []}
	# Gather column for specific patient
	patientColumn = df.loc[:, patientID]
	# Create list of indices to check rows later, instantiate counter
	indicesToProcess, i = [], 0
	# Index last location of common column headers
	end = list(df.columns).index('FORMAT')
	# Check each variant row of the patient data
	for phasingPattern in patientColumn:
		if phasingPattern != '0|0':
			# Parse phasing data, returns '#', '#' [strings]
			leftChromosome, rightChromosome = processPhasedData(phasingPattern)
			# Gather all variant information (i.e. chr#, pos, rsid, etc.)
			variant_df = df.iloc[i, 0:end]
			variant_df = variant_df.to_dict()
			# Sort variants according to L or R from L|R
			if leftChromosome != '0': phasedVariants['ChrL'].append(variant_df)
			if rightChromosome != '0': phasedVariants['ChrR'].append(variant_df)
		# Counter
		i += 1
	return phasedVariants, patientID, df

def buildCSVToQuery(patientIDs, dataframe, csv):
	phasedList = []
	# Create a list of processed phased variants for each patient
	proteinDict = pd.read_csv(csv).to_dict(orient='list')
	chrStrings, proteinStrings = proteinDict['queriedString'], proteinDict['protein_string']
	proteinChanges = dict()
	for i in range(len(chrStrings)):
		protein = proteinStrings[i]
		proteinChanges[chrStrings[i]] = protein

	for patient in patientIDs:
		phasedDict = dict()
		phasedDict['ID'] = patient
		ChrL, ChrR = gatherPhasedVariants(patient, dataframe)[0]['ChrL'], \
		gatherPhasedVariants(patient, dataframe)[0]['ChrR']
		# Add chromosome info before manipulation
		phasedDict['ChrL'], phasedDict['ChrR'] = ChrL, ChrR
		ChrLCopy = copy.deepcopy(ChrL)
		ChrRCopy = copy.deepcopy(ChrR)
		ChrOps = [ChrLCopy, ChrRCopy]
		for i in range(len(ChrOps)):
			proteinOps = ['LeftProtein', 'RightProtein']
			newData = []
			for variant in ChrOps[i]:
				chrFormatted = formatChrAndPos(variant['#CHROM'], variant['POS'], \
					variant['REF'], variant['ALT'])
				variantProtein = proteinChanges[chrFormatted]
				if type(variantProtein) == float:
					variantProtein = None
				newData.append(variantProtein)
			phasedDict[proteinOps[i]] = newData
		phasedList.append(phasedDict)
	# Create a dataframe to hold processed information
	# TODO: Add querying, smoking/alcohol use
	data = phasedList
	df = pd.DataFrame(data, columns = ['ID', 'ChrL', 'ChrR', 'LeftProtein', 'RightProtein'])
	# write dataframe to CSV for human readable format
	df.to_csv('outputs/processed.csv', index = False)
	return df

def grabVariant(chrNumAndPosition):
	# Example input: 'chr7:g.117589482A>G'
	geneDict, overallDict, protein_result = None, None, None
	# Build myvariant object
	mv = myvariant.MyVariantInfo()
	# Query to select dbsnp dictionary from overall dictionary
	if mv:
		try:
			dbsnpDict = mv.getvariant(chrNumAndPosition, assembly='hg38')['dbsnp']
			overallDict = mv.getvariant(chrNumAndPosition, assembly='hg38')
			if dbsnpDict:
				# Gather gene information from dbSNP dictionary
				geneDict = dbsnpDict['gene']
				protein_result = findCommonProteinName(overallDict)
				print('Search Query Successful')
		except TypeError:
			print('Search Query Not Found')
	return geneDict, overallDict, protein_result

def formatChrAndPos(chromNum, pos, ref, alt):
	# Inputs: strings from columns in Beagle dataset
	formattedString = ('chr{}:g.{}{}>{}').format(chromNum, pos, ref, alt)
	# Output string ready to be processed in grabVariant()
	return formattedString

def addQueryStrings(processedDf):
	for groupOfVariants in processedDf['ChrL']:
		for variant in groupOfVariants:
			chromNum = str(variant['#CHROM'])
			pos = str(variant['POS'])
			ref = str(variant['REF'])
			alt = str(variant['ALT'])
			queryString = formatChrAndPos(chromNum, pos, ref, alt)
	return 

def addDBSNPInfo(df, path):
	# Error in that it returns 'V470M' instead of 'M470V'
	# Build List of Queried Information from dbSNP
	queriedList = []
	dbSNPList = []
	generalList = []
	proteinList = []
	for i in range(len(df['#CHROM'])):
		chromNum = df['#CHROM'][i]
		pos = df['POS'][i]
		ref = df['REF'][i]
		alt = df['ALT'][i]
		toQuery = formatChrAndPos(chromNum,pos,ref,alt)
		dbsnpVar, genVar, protein = grabVariant(toQuery)[0], grabVariant(toQuery)[1], grabVariant(toQuery)[2]
		queriedList.append(toQuery)
		dbSNPList.append(dbsnpVar)
		generalList.append(genVar)
		if protein:
			protein = protein.split('|')[1].strip().replace('p.','')
		proteinList.append(protein)

	df2 = pd.DataFrame()
	df2['queriedString'] = queriedList
	#df2['dbsnp'] = dbSNPList
	#df2['overall'] = generalList
	df2['protein_string'] =  proteinList
	df2.to_csv(path, index = False)
	return df2

def findCommonProteinName(overallDict):
	result = None
	try:
		result = overallDict['emv']['egl_protein']
		print('Protein Query Successful')
	except KeyError:
		print('Protein Query Not Found')
	return result

def compareQueriedVariantsToVariantDatabase(queried_csv, database_csv):
	# Input: csv file with variants in database and queried list
	# Output: csv file with 'Present in Database' and 'Not Present in Database'
	dfQueried = pd.read_csv(queried_csv)
	dfDatabase = pd.read_csv(database_csv)
	# Prepare lists to read
	unprocessedQueriedList = dfQueried['protein_string'].tolist()
	unprocessedDatabaseList = dfDatabase['Variant'].tolist()
	databaseList, queriedList = [], []
	# Remove NaN values
	for item in unprocessedDatabaseList:
		if type(item) == str:
			databaseList.append(item)
	for item in unprocessedQueriedList:
		if type(item) == str:
			queriedList.append(item)
	# Create sets to compare
	queriedSet = set(queriedList)
	databaseSet = set(databaseList)
	print('queried', queriedSet)
	# Compare Sets
	presentSet = queriedSet & databaseSet
	notPresentSet = queriedSet - presentSet
	print('present', presentSet)
	print('not present', notPresentSet)
	return

def buildPatientDataframeToGraph(processed_csv, all_variants_csv, CRF_csv):
	# Note bug where G551D must be added in manually b/c not recognized
	# by automated annotation script

	# Read CSV of phased variant data
	processed_df = pd.read_csv(processed_csv)
	# Read current database of variant information
	variant_info = pd.read_csv(all_variants_csv, index_col = 0)
	# Read CRF data
	crf_df = pd.read_csv(CRF_csv)

	right = processed_df['RightProtein']
	left = processed_df['LeftProtein']
	ids = processed_df['ID']

	cftr_df = generateEffectiveCFTRFunction(left, right, ids, variant_info)
	cftr_df.set_index('ID')

	return crf_df, cftr_df

def generateEffectiveCFTRFunction(chrL, chrR, ids, variant_info):
	variant_set = set(variant_info.index.tolist())
	output = []
	for i in range(len(chrR)):
		# String Processing
		left_options = chrL[i].replace('[','').replace(']','').replace("'", '').replace(" ", '').split(',')
		right_options = chrR[i].replace('[','').replace(']','').replace("'", '').replace(" ", '').split(',')
		## To do establish overall impact
		overall_impact_LR = {'Residual':[],'Ivocaftor':[],'Lumacaftor':[],'Combination Therapy':[]}
		for chrList in [left_options, right_options]:
			var_dict = {'Residual':None,'Ivocaftor':None,'Lumacaftor':None,'Combination Therapy':None}
			for key in var_dict.keys():
				variant_impact_list = []
				for item in chrList:
					# Only check variants that quantitative data is available for
					individual_impact = 100
					if item in variant_set:
						fxn = variant_info.loc[item][key]
						# if evidence of quantitative impact add to list
						if fxn != 'None':
							variant_impact_list.append(float(fxn))
						# else assume WT function
						else:
							variant_impact_list.append(100)
					else:
						variant_impact_list.append(100)
				# Add one chromosome impact
				if variant_impact_list != []:
					overall_impact_LR[key].append(statistics.mean(variant_impact_list)/2)
		# Take sum each half of the contribution to overall function
		for key in overall_impact_LR.keys():
			overall_impact_LR[key] = sum(overall_impact_LR[key])
		# Add ID key
		overall_impact_LR['ID'] = ids[i]
		output.append(overall_impact_LR)
	outputdf = pd.DataFrame.from_dict(output)
	outputdf = outputdf[['ID', 'Residual', 'Ivocaftor', 'Lumacaftor', 'Combination Therapy']]
	return outputdf

def generate_graphs_for_df_of_patients(crf_df, cftr_df, patientgroup):
	# Send all outputs to the same subdirectory
	path_prefix = 'outputs/'+ patientgroup
	# Base Case for one patient
	for i in range(len(crf_df['ID'])):
		crf_data = crf_df.iloc[i]
		fxn_data = cftr_df.iloc[i]
		print(crf_data, fxn_data)

		# Check that Patient ID's match from separate CSVs
		if crf_data['ID'] == fxn_data['ID']:
			# Name file and graph
			name = crf_data['ID'].replace('.pjt','')
			pt_subfolder = path_prefix + '/' + name + '/'
			for item in ['Residual', 'Ivocaftor', 'Lumacaftor', 'Combination Therapy']:
				# Establish path extension
				path_extension = name + '_output_' + item
				# Assign datapoint
				datapoint = fxn_data[item]
				# Check if data point is not 100 (no point in graphing WT)
				if not math.isclose(100, datapoint, abs_tol=1e-8):
					# Create a folder for patient if not present
					if not os.path.exists(pt_subfolder):
						os.makedirs(pt_subfolder)
					# CFTR Function
					effective_cftr_fxn = datapoint
					# Initiate Mechanistic Model
					cell_sim = Duct_Cell()
					# Add Variant Impact
					cell_sim.add_variant_influence(effective_cftr_fxn/100)
					# Assign env. vars
					smoking_ever = crf_data['Smoking']
					smoking_status = crf_data['Current Smoker / Past Smoker']
					alc_ever = crf_data['Alcohol']
					alcohol_status = crf_data['Current Drinker / Past Drinker']
					# Patient Status at Baseline
					if item == 'Residual':
						# Add Smoking / Drinking
						if alc_ever == "TRUE":
							cell_sim.add_alcohol_influence(0.8)
						if smoking_ever == "TRUE":
							if smoking_status == "Current":
								cell_sim.add_smoking_influence(5/11)
							elif smoking_status == "Past":
								cell_sim.add_smoking_influence(8.5/11)
					

					# Run model
					
					# Instantiate New Cell for Patient
					wt_cell = Duct_Cell()
					input_data = wt_cell.input_dict

					init = copy.deepcopy(init_cond)

					# Construct Arrays from Duct Model System Equations
					def generate_source_array(input_data, key):
						choice_array = run_model_CFTR(input_data, 20000, 120000, 200000)[0][key]
						if key == 'time':
							choice_array /=  20000
						return choice_array

					# Individual graphing arrays for ColumnDataSource
					wt_bi_l = generate_source_array(init, 'bl')
					wt_bi_i = generate_source_array(init, 'bi')
					pt_bi_l = generate_source_array(cell_sim.input_dict, 'bl')
					pt_bi_i = generate_source_array(cell_sim.input_dict, 'bi')
					wt_cl_l = 160 - wt_bi_l
					wt_cl_i = generate_source_array(init, 'ci')
					pt_cl_l = 160 - pt_bi_l
					pt_cl_i = generate_source_array(cell_sim.input_dict, 'ci')
					time = generate_source_array(input_data, 'time')

					item_title = item
					if item == 'Ivocaftor':
						item_title = 'Ivacaftor'
					title =  'Patient Information for ' + name + ' (' + item_title +') '
					if item == 'Residual':
						if alc_ever == "TRUE"  and smoking_ever	== "TRUE":
							title += '- Alcohol and Tobacco Influences Added'
						elif alc_ever == "TRUE" and smoking_ever != "TRUE":
							title += '- Alcohol Influence Added'
						elif alc_ever != "TRUE" and smoking_ever == "TRUE":
							title += '- Tobacco Influence Added'
						elif alc_ever != "TRUE" and smoking_ever != "TRUE":
							title += '- No Tobacco or Alcohol Use Reported in CRF'
					else:
						title += '- In addition to abstaining from Tobacco & Alcohol'

					# Bicarb Plot
					bi_legend = Legend(items=[])
					bi_plot = figure(plot_height=800, plot_width=1000,  title=title)
					bi_plot.add_layout(bi_legend, 'right')
					bi_plot.legend.click_policy = 'hide'
					bi_plot.line(time, pt_bi_l, line_width=3, line_color = '#34344A', legend="Patient Luminal HCO3-")
					bi_plot.line(time, pt_bi_i, line_width=3, line_color = '#7FE0CB', legend="Patient Intra HCO3-")
					bi_plot.line(time, wt_bi_l, line_width=3, line_color = '#34344A', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Luminal HCO3-")
					bi_plot.line(time, wt_bi_i, line_width=3, line_color = '#7FE0CB', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Intracellular HCO3-")
					bi_plot.legend.label_text_font = 'gilroy'
					bi_plot.title.text_font = 'gilroy'
					bi_plot.title.text_font_style = 'bold'
					bi_plot.yaxis.axis_label_text_font = 'gilroy'
					bi_plot.yaxis.axis_label_text_font_style = 'normal'
					bi_plot.xaxis.axis_label_text_font = 'gilroy'
					bi_plot.xaxis.axis_label_text_font_style = 'normal'

					bi_plot.yaxis.axis_label = 'Bicarb Conc. (mM)'
					bi_plot.xaxis.axis_label = 'Time'

					# Chloride Plot
					cl_legend = Legend(items=[])
					cl_plot = figure(plot_height=800, plot_width=1000,  title=title)
					cl_plot.add_layout(cl_legend, 'right')
					cl_plot.legend.click_policy = 'hide'
					cl_plot.line(time, pt_cl_l, line_width=3, line_color = '#34344A', legend="Patient Luminal Cl-")
					cl_plot.line(time, pt_cl_i, line_width=3, line_color = '#7FE0CB', legend="Patient Intra Cl-")
					cl_plot.line(time, wt_cl_l, line_width=3, line_color = '#34344A', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Luminal Cl-")
					cl_plot.line(time, wt_cl_i, line_width=3, line_color = '#7FE0CB', alpha=0.25, line_dash='dashed', line_cap='round', legend="WT Intracellular Cl-")
					cl_plot.legend.label_text_font = 'gilroy'
					cl_plot.legend.background_fill_color = '#F4F1E1'
					cl_plot.legend.background_fill_alpha = 0.25
					cl_plot.title.text_font = 'gilroy'
					cl_plot.title.text_font_style = 'bold'
					cl_plot.yaxis.axis_label_text_font = 'gilroy'
					cl_plot.yaxis.axis_label_text_font_style = 'normal'
					cl_plot.xaxis.axis_label_text_font = 'gilroy'
					cl_plot.xaxis.axis_label_text_font_style = 'normal'

					cl_plot.yaxis.axis_label = 'Chloride Conc. (mM)'
					cl_plot.xaxis.axis_label = 'Time'

					save(bi_plot, filename= pt_subfolder + path_extension + '_hco3' + '.html')
					save(cl_plot, filename= pt_subfolder + path_extension + '_cl' + '.html')


		else:
			print('Patient ID\'s do not match')

	return



# Perform functions
# df1 = appendColumnHeaders('colnames.txt', 'beagle_out_CFTR_genotype.txt')
# #addDBSNPInfo(df1, 'outputs/queriedVariantList.csv')
# exportExcel('outputs/beagle_joined.xlsx', df1)
# patientIDs = gatherPatientIDs(df1)
# #df2 = buildCSVToQuery(patientIDs, df1, 'outputs/queriedVariantList.csv')


crf_df, cftr_df = buildPatientDataframeToGraph('outputs/processed.csv', 'inputs/all_variants.csv', 'inputs/patient_env_choices.csv')
generate_graphs_for_df_of_patients(crf_df, cftr_df, 'hopkins_patients')










