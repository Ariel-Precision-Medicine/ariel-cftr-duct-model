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

def buildCSVToQuery(patientIDs, dataframe):
	phasedList = []
	# Create a list of processed phased variants for each patient
	for patient in patientIDs:
		phasedDict = dict()
		phasedDict['ID'] = patient
		ChrL, ChrR = gatherPhasedVariants(patient, dataframe)[0]['ChrL'], \
		gatherPhasedVariants(patient, dataframe)[0]['ChrR']
		phasedDict['ChrL'], phasedDict['ChrR'] = ChrL, ChrR
		phasedList.append(phasedDict)
	# Create a dataframe to hold processed information
	# TODO: Add querying, smoking/alcohol use
	data = phasedList
	df = pd.DataFrame(data, columns = ['ID', 'ChrL', 'ChrR'])
	# write dataframe to CSV for human readable format
	df.to_csv('outputs/processed.csv')
	return df

def grabVariant(chrNumAndPosition):
	# Example input: 'chr7:g.117589482A>G'
	geneDict, overallDict = None, None
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
				print('Search Query Successful')
		except TypeError:
			print('Search Query Not Found')
	return geneDict, overallDict

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

def addDBSNPInfo(df):
	# Build List of Queried Information from dbSNP
	queriedList = []
	dbSNPList = []
	generalList = []
	for i in range(len(df['#CHROM'])):
		chromNum = df['#CHROM'][i]
		pos = df['POS'][i]
		ref = df['REF'][i]
		alt = df['ALT'][i]
		toQuery = formatChrAndPos(chromNum,pos,ref,alt)
		dbsnpVar, genVar = grabVariant(toQuery)[0], grabVariant(toQuery)[1]
		queriedList.append(toQuery)
		dbSNPList.append(dbsnpVar)
		generalList.append(genVar)
	df2 = pd.DataFrame()
	df2['queriedString'] = queriedList
	df2['dbsnp'] = dbSNPList
	df2['overall'] = generalList
	df2.to_csv('outputs/queriedVariantList.csv')
	return df2

# Perform functions
df1 = appendColumnHeaders('colnames.txt', 'beagle_out_CFTR_genotype.txt')
addDBSNPInfo(df1)
exportExcel('outputs/beagle_joined.xlsx', df1)
patientIDs = gatherPatientIDs(df1)
df2 = buildCSVToQuery(patientIDs, df1)

#grabVariant('chr7:g.117509093G>A')
#addQueryStrings(df2)


