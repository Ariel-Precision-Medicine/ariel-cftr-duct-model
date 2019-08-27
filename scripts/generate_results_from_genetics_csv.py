# generate_results_from_genetics_csv.py

# Cameron Breze for Ariel Precision Medicine

# Purpose of this script is to generate sample data to test
# automation of report generation for basic UPMC dataset and
# future work

from oop_duct_model import Duct_Cell
import pandas, random, secrets

# Initialize Cell
sampleCell = Duct_Cell()

# Create List of Possible Variants
cftrVariantOptions = sampleCell.gen_var_menu()

# Create random list of variants of a specific length for test case
def generateRandomList(desiredLength):
	lengthOptions = len(cftrVariantOptions)
	doubleLengthOptions = 2 * lengthOptions
	cftrWeights = [1/(doubleLengthOptions)] * lengthOptions
	cftrWeights[0] = cftrWeights[0] * (lengthOptions + 1)
	return random.choices(cftrVariantOptions, 
		weights = cftrWeights, k = desiredLength)

# Create random list of pt. identifiers of a specific length for test case
def generateRandomPatientList(desiredLength):
	patientList = []
	for i in range(desiredLength):
		patientList.append(secrets.token_urlsafe(12))
	return patientList

# Generate dataframe of random sample variant in patient population
def generateRandomDataFrame(desiredLength):
	patientList = generateRandomPatientList(desiredLength)
	varList1 = generateRandomList(desiredLength)
	varList2 = generateRandomList(desiredLength)
	data = {'Patient ID': patientList, 'Variant 1':varList1, 'Variant 2':varList2}
	dataFrame = pandas.DataFrame(data)
	return dataFrame

# Generate CSV (Excel sheet) with random patient information, writes CSV
def buildCSV(desiredLength, filename):
	df = generateRandomDataFrame(desiredLength)
	df.to_csv(filename, index = False)
