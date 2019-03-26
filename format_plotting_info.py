#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

#returnMode spex
#slightly a misnomer finds the nth biggest non signature of interest and returns it
#n is 1 indexed for readability so we subtract 1
def find_nth_biggest_non_primary_sig(row, primarySig, n, prefix='Signature.'):
    colNames = row.to_dict().keys()

    signatureColumns = [i for i in colNames if prefix in i]
    rowSigsOnly = row[signatureColumns]
    rowAsDict = rowSigsOnly.to_dict()
    items = rowAsDict.items()
    sortedItems = sorted(items, key=lambda x: x[1], reverse=True)

    #remove the dominant signature
    sortedItemsLimited = [z for z in sortedItems if z[0] != primarySig]
    return sortedItemsLimited[n - 1][0], sortedItemsLimited[n - 1][1]

def merge_signature_columns(df, mode='Stratton', drop=True, smokingMerge=False, confidence=True, mean=True, prefix='mean_'):
	if mode == 'Stratton':
		if confidence: df['confidence_APOBEC'] = df.apply(lambda row: max(row['confidence_2'], row['confidence_13']), axis=1)
		if mean: df[prefix + 'APOBEC'] = df.apply(lambda row: row[prefix + '2'] + row[prefix + '13'], axis=1)
		if confidence: df['confidence_MMR'] = df.apply(lambda row: max(row['confidence_6'], row['confidence_15'], row['confidence_20'], row['confidence_21'], row['confidence_26']), axis=1)
		if mean: df[prefix + 'MMR'] = df.apply(lambda row: row[prefix + '6'] + row[prefix + '15'] + row[prefix + '20'] + row[prefix + '21'] + row[prefix + '26'], axis=1)	
		if mean:
			if smokingMerge: #smoking merge, if specified merges the smoking signature, mutyh, aflatoxin and chewing tobacco (which are all usually smoking)
				df[prefix + 'SMOKING'] = df.apply(lambda row: row[prefix + '4'] + row[prefix + '18'] + row[prefix + '24'] + row[prefix + '29'], axis=1)
		dropCols = []
		if mean: dropCols += [prefix + '2', prefix + '13', prefix + '6', prefix + '15', prefix + '20', prefix + '21', prefix + '26']
		if smokingMerge: dropCols += [prefix + '4', prefix + '18', prefix + '24', prefix + '29']
		if confidence: dropCols += ['confidence_2', 'confidence_13', 'confidence_6', 'confidence_15', 'confidence_20', 'confidence_21', 'confidence_26']
		if drop: #drop cols if asked
			df = df.drop(dropCols, axis=1)
		return df 

#returns the dominant signautre for a row of a df expressed as a dict with the specified signatures under consideration
def get_dominant_signature(rowAsDict, prefix='mean'):
	cols = [x for x in rowAsDict.keys() if 'Signature.' in x]
	tupList = []
	for key, value in rowAsDict.items():
		if prefix in key: tupList.append((key, value))
	sortedSigs = sorted(tupList, key = lambda tup: tup[1], reverse=True)
	return sortedSigs[0][0]

def calculate_ordering_val(row):
	return 0


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--signature_data_file', help='file with signature information', default='/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
	parser.add_argument('--writeFileName', help='filename for file', default='signatureWaterfallPlot')
	parser.add_argument('--writeFileDir', help='directory to write file in', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles')
	parser.add_argument('--idsFile', help='path to file listing ids to analyze', default='test_signature_case_ids.txt')
	parser.add_argument('--primarySortSig', help='The signature to order/display data by', default='Signature.AGE')
	parser.add_argument('--minMutThreshold', help='The minumum number of mutations below which we wont call signatures', default=10)

	args = parser.parse_args()

	#READ in different formatting info

	#ids 
	caseIdsFile = open(args.idsFile)
	lines = caseIdsFile.readlines()
	caseIds = set([i.strip(',').strip('\n').strip('\t') for i in lines]) #supports tsv, csv or newline delineated files

	#ordering of signatures information
	sigOrderingInfo = open('signatureOrdering.txt').readlines()[0].strip('\n')
	sigOrderingDict =  dict(zip([i for i in sigOrderingInfo.split(',')], [i for i in range(len(sigOrderingInfo.split(',')))]))

	#load in the signatures information and add some utility columns
	print 'loading in signatures information'
	signaturesDf = pd.read_table(args.signature_data_file)
	signaturesDf = signaturesDf[signaturesDf['Tumor_Sample_Barcode'].isin(caseIds)] #only include the ids we want
	signaturesDf['pid'] = signaturesDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	signaturesDf['tid'] = signaturesDf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])

	#see if the prefix is mean or signature
	#rename columns as we need
	signatureRenameDict = dict(zip(['mean_' + str(i) for i in range(1,31)], ['Signature.' + str(i) for i in range(1,31)]))
	signaturesDf = signaturesDf.rename(columns= signatureRenameDict)
	signaturesDf = merge_signature_columns(signaturesDf, mode='Stratton', prefix='Signature.')
	signaturesDf = signaturesDf.rename(columns={'Signature.1': 'Signature.AGE',
												'Signature.3': 'Signature.HRD',
												'Signature.4': 'Signature.SMOKING',
												'Signature.7': 'Signature.UV',
												'Signature.10': 'Signature.POLE',
												'Signature.11': 'Signature.TMZ',
												'Signature.14': 'Signature.POLE_plus_MMR'}) #rename a second time to more human readable signature names


	#create a column for the dominant signatures
	signaturesDf['signatureOfInterestMagnitude'] = signaturesDf[args.primarySortSig]
	signaturesDf['signatureOfInterestName'] = args.primarySortSig
	signaturesDf['dominantSignature'] = signaturesDf.apply(lambda row: get_dominant_signature(row.to_dict(), prefix='Signature.'), axis=1) #get the dominant signature

	#get the second most common signature
	signaturesDf['secondPredominantSigName'] = signaturesDf.apply(lambda row: find_nth_biggest_non_primary_sig(row, args.primarySortSig, n=1, prefix='Signature.')[0], axis=1)
	signaturesDf['secondPredominantSigMagnitude'] = signaturesDf.apply(lambda row: find_nth_biggest_non_primary_sig(row, args.primarySortSig, n=1, prefix='Signature.')[1], axis=1)
	signaturesDf['secondPredominantSigName'] = signaturesDf['secondPredominantSigName'].apply(lambda x: x if x in sigOrderingDict.keys() else 'other')

	#get the third most common signature
	signaturesDf['thirdPredominantSigName'] = signaturesDf.apply(lambda row: find_nth_biggest_non_primary_sig(row, args.primarySortSig, n=2, prefix='Signature.')[0], axis=1)
	signaturesDf['thirdPredominantSigMagnitude'] = signaturesDf.apply(lambda row: find_nth_biggest_non_primary_sig(row, args.primarySortSig, n=2, prefix='Signature.')[1], axis=1)
	signaturesDf['thirdPredominantSigName'] = signaturesDf['thirdPredominantSigName'].apply(lambda x: x if x in sigOrderingDict.keys() else 'other')

	#creates an "ordering value which R will use to line up the plot---lower means closer to the left of the screen"
	signaturesDf['orderingVal'] = signaturesDf.apply(lambda row:
		0 - row['signatureOfInterestMagnitude'] if row['dominantSignature'] == row['signatureOfInterestName']
		else sigOrderingDict[row['dominantSignature']] - row['secondPredominantSigMagnitude'] if row['dominantSignature'] in sigOrderingDict
		else 100 - row['secondPredominantSigMagnitude']
		, axis=1)
		#sigOrderingDict[row['dominantSignature']] if row['dominantSignature'] in sigOrderingDict else len(sigOrderingDict.keys()) + 1, axis = 1)

	#FORMAT TO IGNORE the cases with not enough mutations
	signaturesDf['orderingVal'] = signaturesDf.apply(lambda row: row['orderingVal'] if row['Nmut'] >= args.minMutThreshold else 1000 - row['Nmut'], axis=1)
	signaturesDf['signatureOfInterestMagnitude'] = signaturesDf.apply(lambda row: row['signatureOfInterestMagnitude'] if row['Nmut'] >= args.minMutThreshold else None, axis=1)
	signaturesDf['secondPredominantSigMagnitude'] = signaturesDf.apply(lambda row: row['secondPredominantSigMagnitude'] if row['Nmut'] >= args.minMutThreshold else None, axis=1)
	signaturesDf['thirdPredominantSigMagnitude'] = signaturesDf.apply(lambda row: row['thirdPredominantSigMagnitude'] if row['Nmut'] >= args.minMutThreshold else None, axis=1)

	writeFilePath = os.path.join(args.writeFileDir, args.writeFileName + '.csv')
	print 'writing file to ', writeFilePath
	signaturesDf.to_csv(writeFilePath, sep='\t', index=False)
	

if __name__ == '__main__':
    main()















