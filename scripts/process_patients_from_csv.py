# process_patients_from_csv.py

# Cameron Breze for Ariel Precision Medicine

# Purpose of this script is to process patient CFTR variants
# and return a pandas Dataframe with:

''' 
| Patient ID | - string
| % WT Function (overall) | - float
| % Change in peak ion flux | - float
| VarA | - dict of common name, HGVS, cDNA, 
	likelihood of drug response & % function w.r.t WT
| VarB | - dict of common name, HGVS, cDNA, 
	likelihood of drug response & % function w.r.t WT
| bokeh | - dict of graph objects containing HCO3- and Cl- transport info
'''

import pandas

# String processing (removes all characters after plus/minus)
def remove_standard_error(text_string):
		return text_string[0:text_string.find(u'\u00b1')].strip()

# Create variant dictionary with relevant information
def buildVariantDict(variantString, csv):
	# Build dictionary
	varDict = dict()
	# Create master df from csv
	df = pandas.read_csv(csv)
	# Create single row sub-df to collect individual data points
	sub_df = df.loc[df.Variant == variantString]
	# Collect variant information and match to dictionary keys
	varDict['Name'] = variantString
	varDict['HGVS'] = sub_df.iloc[0]['HGVS']
	varDict['cDNA'] = sub_df.iloc[0]['cDNA']
	varDict['% WT function'] = \
		float(remove_standard_error(sub_df.iloc[0]['Residual']))
	varDict['Lumacaftor'] = \
		float(remove_standard_error(sub_df.iloc[0]['Lumacaftor']))
	varDict['Ivocaftor'] = \
		float(remove_standard_error(sub_df.iloc[0]['Ivocaftor']))
	varDict['Ivocaftor and Lumacaftor'] = \
		float(remove_standard_error(sub_df.iloc[0]['Ivocaftor and Lumacaftor']))
	return varDict



print(buildVariantDict('F311L', 'cutting_variant_data.csv'))
print(buildVariantDict('I336K', 'cutting_variant_data.csv'))