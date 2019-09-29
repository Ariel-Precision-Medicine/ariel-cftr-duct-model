# beagle_text_processing.py

'''
Cameron Breze for Ariel Precision Medicine
Purpose: Processed Phased CFTR data from Beagle Output

Dependencies:
-pandas
-openpyxl

'''

import pandas as pd

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
			# Sort variants according to L or R from L|R
			if leftChromosome != '0': phasedVariants['ChrL'].append(variant_df)
			if rightChromosome != '0': phasedVariants['ChrR'].append(variant_df)
		# Counter
		i += 1
	return phasedVariants, patientID, df

def buildCSVToQuery(patientIDs, dataframe):
	phasedList = []
	# Create a list of processed phased variants for each patient
	for patient in patientIDs:
		phasedList.append(gatherPhasedVariants(patient, dataframe)[0])
	# Create a dataframe to hold processed information
	# TODO: Add querying, smoking/alcohol use
	data = {'patientIDs': patientIDs, 'phasedVariants': phasedList}
	df = pd.DataFrame(data)
	# write dataframe to CSV for human readable format
	df.to_csv('outputs/processed.csv')
	return df


# Perform functions
df = appendColumnHeaders('colnames.txt', 'beagle_out_CFTR_genotype.txt')
exportExcel('outputs/beagle_joined.xlsx', df)
patientIDs = gatherPatientIDs(df)
buildCSVToQuery(patientIDs, df)


