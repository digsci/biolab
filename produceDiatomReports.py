import argparse
import pandas
import re
from collections import defaultdict

def main():
	options = parseArguments()
	resultsFile = options.folder + "/visualised_taxonomy/otu_table.diatomsonly_L1.txt"
	countsFile = options.folder + "/diatomSequenceCounts.txt"

	# Import the sequence counts
	counts = {}
	for line in open(countsFile, "rU"):
		linelist = re.split(' ',line)
		filename = linelist[0]
		count = int(45367) #int(linelist[1].rstrip())
		#filename will be sample.passedQC.fastq
		filenamesplit = re.split('\.',filename)
		sample = filenamesplit[0].replace(options.folder,'').replace('/','')
		counts[sample] = count

	# Import the lookup if there
	if (options.lookup):
		lookup = {}
		# FERA sample name \t EA sample name
		for line in open(options.lookup, "rU"):
			linelist = re.split('\t',line)
			lookup[linelist[0]] = linelist[1].rstrip()

	# Import the results into a pandas dataframe
	data = pandas.read_csv(resultsFile,header=1,sep='\t',index_col=0)
	data.index.names = ['Taxonomy']

	# screen on sequence counts. If a sequence count is <= 3000 then the sequencing will need repeating
	# add rag row to the dataframe
	repeats = []
	passes = []
	for sample in counts:
		count = int(counts[sample])
		if count <= 3000:
			repeats.append(sample)
			if count > 3000:
				passes.append(sample)
	passed = data.loc[passes]
	repeat = data.loc[repeats]

	# change the dataframe column names (sample names from lookup) if required.
	# if a name isn't in the lookup then it won't change it and keep the fera one.
	passed.rename(index=str, columns=lookup, copy=True, inplace=True)
	repeat.rename(index=str, columns=lookup, copy=True, inplace=True)
	# Export the dataframe back out to csv
	passed.to_csv(options.folder + "/Abundances.pass.csv")
	repeat.to_csv(options.folder + "/Abundances.fail.csv")
        print "Reports completed"

def parseArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('--folder', help='The folder containing the inputs for diatomPipeline.sh and the files produced by the pipeline', required=True)
	parser.add_argument('--lookup', help='Optional. If there is a lookup of samples names to EA barcodes. Tab-delimited file where column 1 is the FERA reference as described in the MiSeq sample sheet; and column 2 is the EA barcode number')
	return parser.parse_args()

if __name__ == '__main__':
	main()
