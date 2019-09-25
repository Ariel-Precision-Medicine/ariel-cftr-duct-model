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

def gatherPhasedVariants(patientID, df):
	# Gather column for specific patient
	patientColumn = df.loc[:, patientID]
	# Create variant subset with none result removed
	indicesToProcess = []
	for phasingPattern in patientColumn:
		if phasingPattern == '0|0':
			print('None')
		else:
			print('Something')
			
	print(patientColumn)
	return





# Perform functions
df = appendColumnHeaders('colnames.txt', 'beagle_out_CFTR_genotype.txt')
exportExcel('outputs/beagle_joined.xlsx', df)
gatherPatientIDs(df)

gatherPhasedVariants('GS_AR18Q10009_V4.pjt', df)

