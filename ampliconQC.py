#!/usr/bin/python

from __future__ import division
import argparse
import glob
import re
import os
import subprocess
import pygal
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
from Bio import SeqIO
from collections import Counter

def main():
	options = processArguments()

	# Get the real path of the directory. User may have given a relative path.
	dataDirectory = os.path.realpath(options.data) + "/"

	# Rename the files to the format sampleName.readDirection.fastq.gz
	renamed = renameFiles(dataDirectory)

	# Move the 'undetermined' files to the raw_data folder as
	moved = moveUndetermined(dataDirectory)

	# Trim off the primers from the sequences
	trimmed = runCutadapt(options,dataDirectory)

	# Run sickle paired-end to trim off bad quality 3' bases
	peQC = runSicklePE(options,dataDirectory)

	# Run pear to merge the R1 and R2 reads
	joined = runPear(options,dataDirectory)

	# Run sickle single-end to remove post-merging bad quality sequences
	seQC = runSickleSE(options,dataDirectory)

	# If QIIME prep is required:
	if (options.qiime == True):
		prepped = runQiimePrep(options,dataDirectory)

	# If histograms are requested:
	if (options.histograms == True):
		sequenceFiles = glob.glob(dataDirectory + "*.passedQC.fastq")
		for sequenceFile in sequenceFiles:
			graph = drawHistogram(sequenceFile)

	# FINISHED
	print "FINISHED Amplicon quality control!"

# Draw a histogram given a file of sequences
def drawHistogram(sequenceFile):
	filename, extension = os.path.splitext(sequenceFile)
	extension = re.sub('\.','',extension) #replace the full stop in the ext.
	outfile = str(sequenceFile) + ".histogram.svg"
	sequences = [len(rec) for rec in SeqIO.parse(sequenceFile,extension)]
	counts = Counter(sequences)
	plot = pygal.XY(show_x_guides=True,show_legend=False,title="Sequence Length Histogram: "+sequenceFile, x_title="Sequence length (nt)", y_title="Number of sequences")
	plot.add('Sequence Lengths',counts.items(),show_dots=False)
	plot.render_to_file(outfile)
	print "Histogram " + outfile + " finished"
	return True

# This function changes the sample name to a QIIME-acceptable one by firstly changing any '-' in the
# sample name to '.'. It then adds the modified sample name to each sequence in the sample fastq file.
# All the modified sample files are then concatenated into one large file and an incremental number applied to each
# sequence in the file. The temporary files are removed, and a skeleton mapping file is produced from the modified sample names.
def runQiimePrep(options,directory):
	# Retrieve a list of all the passedQC files
	sequenceFiles = glob.glob(directory + "*.passedQC.fastq")

	# List to hold the amended sample names to go in the skeleton mapping file
	sampleNames = []

	# Amend the sequence ID of each sequence in each file
	for fastqFile in sequenceFiles:
		filenameList = re.split('\.',fastqFile)
		originalSampleName = filenameList[0].replace(directory,'')
		newSampleName = originalSampleName.replace('-','.')
		sampleNames.append(newSampleName)
		tempOutfile = open(directory + "temp." + newSampleName + ".fastq","w")
		# Amend the sequence id's
		for seq in SeqIO.parse(fastqFile,"fastq"):
			name = seq.id
			seq.id = newSampleName
			SeqIO.write(seq,tempOutfile,"fastq")
		tempOutfile.close()
		# Rename the original file to the new sample name file (so no confusion later with the mapping file)
		rename = os.rename(fastqFile,directory + newSampleName + ".passedQC.fastq")

	# Concatenate all the 'temp.*' files
	amendedSequenceFiles = glob.glob(directory + "temp.*")
	allseqs = open(directory + "allpassedQC.unnumbered.fastq","w")
	for seqFile in amendedSequenceFiles:
		for line in open(seqFile,"rU"):
			allseqs.write(line)
	allseqs.close()

	# Remove all the 'temp.*' files
	for sample in sampleNames:
		tempdel = os.remove(directory + "temp." + sample + ".fastq")

	# Add an incremental number to the allpassedQC.unnumbered.fastq files
	qiimeOutFasta = open(directory + "readyForQiime.allsamples.fasta","w")
	i = 1
	for seq in SeqIO.parse(directory + "allpassedQC.unnumbered.fastq","fastq"):
		name = seq.id
		seq.id = name + "_" + str(i)
		SeqIO.write(seq,qiimeOutFasta,"fasta")
		i = i + 1
	qiimeOutFasta.close()

	# Delete the unnumbered fastq file
	deleted = os.remove(directory + "allpassedQC.unnumbered.fastq")

	# Create a skeleton mapping file
	mappingOut = open(directory + "mappingfile.txt","w")
	mappingOut.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n")
	for sample in sampleNames:
		mappingOut.write(sample + "\t\t\t" + sample + "\n")
	mappingOut.close()

	# If histograms are required, draw one for the entire amplicon qiime file
	if (options.histograms == True):
		graph = drawHistogram(directory + "readyForQiime.allsamples.fasta")

	return True

# Run Sickle SE to QC the joined up reads (tested with v1.33)
def runSickleSE(options,directory):
	minLength = calculateMeanLength(directory) - 50 #fwd primer that was trimmed off
	# Run sickle on each file in the directory (simultanously with parallel)
	sickleSEcommand = "ls " + directory + "*.assembled.fastq | parallel -j " + str(options.threads) + " 'sickle se -f {} -t sanger -o {}.passedQC.fastq -n -l " + str(minLength) + " -q 30'"
	sickleSEchild = subprocess.Popen(str(sickleSEcommand),
									  stdout = subprocess.PIPE,
									  stderr = subprocess.PIPE,
									  universal_newlines = True,
									  shell=(sys.platform!="win32"))
	sickleSEoutput,sickleSEerror = sickleSEchild.communicate()
	print sickleSEoutput,sickleSEerror

	# Remove the assembled pear files and rename the passedQC files
	sequenceFiles = glob.glob(directory + "*.fastq")
	samples = []
	for fastqFile in sequenceFiles:
		filenameList = re.split('\.',fastqFile)
		samples.append(filenameList[0])
	sampleNames = list(set(samples))
	for sampleName in sampleNames:
		rm1 = os.remove(sampleName + ".assembled.fastq")
		rename = os.rename(sampleName + ".assembled.fastq.passedQC.fastq",sampleName + ".passedQC.fastq")
	return True

# Calculate the mean length of amplicons across all PEAR assembled sequence files
def calculateMeanLength(directory):
	mean = None
	means = []
	for fastqfile in glob.glob(directory + "*.assembled.fastq"):
		awk = "awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' " + fastqfile
		awkChild = subprocess.Popen(str(awk),
					stdout = subprocess.PIPE,
					stderr = subprocess.PIPE,
					universal_newlines = True,
					shell=(sys.platform!="win32"))
		awkOutput, awkError = awkChild.communicate()
		output = awkOutput.rstrip()
		thismean = float(output)
		means.append(thismean)
	mean = sum(means) / float(len(means))
	return mean


# Run Pear to join the paired ends (tested with v0.9.10)
def runPear(options,directory):
	# Retrieve a list of unique sample names in the directory
	sequenceFiles = glob.glob(directory + "*.fastq.gz")
	samples = []
	for fastqFile in sequenceFiles:
		filenameList = re.split('\.',fastqFile)
		samples.append(filenameList[0].replace(directory,''))
	sampleNames = list(set(samples))

	# Now run pear on each of the pairs of read files
	for sampleName in sampleNames:
		sampleName = directory + sampleName
		pearCommand = "pear -e 1 -f " + sampleName + ".sickle.trimmed.R1.fastq.gz -r " + sampleName + ".sickle.trimmed.R2.fastq.gz -o " + sampleName + " -j " + str(options.threads)
		pearChild = subprocess.Popen(str(pearCommand),
						 stdout = subprocess.PIPE,
						 stderr = subprocess.PIPE,
						 universal_newlines = True,
						 shell=(sys.platform!="win32"))
		pearOutput,pearError = pearChild.communicate()
		print pearOutput,pearError
		# Remove the pear input files
		rm1 = os.remove(sampleName + ".sickle.trimmed.R1.fastq.gz")
		rm2 = os.remove(sampleName + ".sickle.trimmed.R2.fastq.gz")
		# Remove the unneeded pear output files
		rm3 = os.remove(sampleName + ".discarded.fastq")
		rm4 = os.remove(sampleName + ".unassembled.forward.fastq")
		rm5 = os.remove(sampleName + ".unassembled.reverse.fastq")
	return True

# Move the two 'Undetermined.*.fastq.gz' files to the raw data folder as they're not
# needed but often copied
def moveUndetermined(directory):
	try:
		makeDir = os.mkdir(directory + "raw_data")
	except OSError:
		pass #already exists, carry on.
	try:
		moveR1 = os.rename(directory + "Undetermined.R1.fastq.gz", directory + "raw_data/Undetermined.R1.fastq.gz")
		moveR2 = os.rename(directory + "Undetermined.R2.fastq.gz", directory + "raw_data/Undetermined.R2.fastq.gz")
	except OSError:
		#don't bother doing anything, they're not there.
		return True
	return True

# Run sickle PE (tested with v1.33)
def runSicklePE(options,directory):
	# Retrieve a list of unique sample names in the directory
	sequenceFiles = glob.glob(directory + "*.fastq.gz")
	samples = []
	for fastqFile in sequenceFiles:
		filenameList = re.split('\.',fastqFile)
		samples.append(filenameList[0].replace(directory,''))
	sampleNames = list(set(samples))
	# Now run sickle pe on each of the samples

	for sampleName in sampleNames:
		sicklePEcommand = "sickle pe -f " + directory + sampleName + ".R1.fastq.gz.trimmed.fastq.gz -r " + directory + sampleName + ".R2.fastq.gz.trimmed.fastq.gz -t sanger -o "+directory + sampleName+".sickle.trimmed.R1.fastq.gz -p "+directory + sampleName+".sickle.trimmed.R2.fastq.gz -s "+directory + sampleName+".single.trimmed.fastq.gz -x"
		sicklePEchild = subprocess.Popen(str(sicklePEcommand),
						 stdout = subprocess.PIPE,
						 stderr = subprocess.PIPE,
						 universal_newlines = True,
						 shell=(sys.platform!="win32"))
		sicklePEOutput,sicklePEError = sicklePEchild.communicate()
		print sicklePEOutput,sicklePEError

		# Remove the sickle input files and single trimmed files.
		delr1 = os.remove(directory + sampleName + ".R1.fastq.gz.trimmed.fastq.gz")
		delr2 = os.remove(directory + sampleName + ".R2.fastq.gz.trimmed.fastq.gz")
		delsingle = os.remove(directory + sampleName + ".single.trimmed.fastq.gz")
		# Move the raw data files into their own folder
		move2 = os.rename(directory + sampleName + ".R1.fastq.gz",directory + "raw_data/" + sampleName + ".R1.fastq.gz")
		move3 = os.rename(directory + sampleName + ".R2.fastq.gz",directory + "raw_data/" + sampleName + ".R2.fastq.gz")
	return True

# Run cutadapt (tested with v 1.9.1). The error rate is calculated based on the number of
# degeneracies present in the primer sequences (cutadapte treats degeneracies as errors)
def runCutadapt(options, directory):
	errorRate = calculateCutadaptError(options.forward,options.reverse)
	# Run cutadapt on each file in the directory (simultaneously with parallel)
	cutadaptCommand = "ls " + directory + "*.fastq.gz | parallel -j " + str(options.threads) + " 'cutadapt -e " + str(errorRate) + " -b " + str(options.forward) + " -b " + str(options.reverse) + " -o {}.trimmed.fastq.gz {}'"
	cutadaptChild = subprocess.Popen(str(cutadaptCommand),
					 stdout = subprocess.PIPE,
					 stderr = subprocess.PIPE,
					 universal_newlines = True,
					 shell=(sys.platform!="win32"))
	cutadaptOutput,cutadaptError = cutadaptChild.communicate()
	print cutadaptOutput,cutadaptError
	return True


# Calculate the cutadapt error rate to use for this primer pair. Whichever primer has the
# highest error rate is the rate to use (default is 0.1)
def calculateCutadaptError(forward,reverse):
	errorRate = 0.1
	forward = forward.upper()
	reverse = reverse.upper()
	fwdBases = forward.count('A') + forward.count('C') + forward.count('G') + forward.count('T')
	revBases = reverse.count('A') + reverse.count('C') + reverse.count('G') + reverse.count('T')
	fwdDegeneracies = len(forward) - fwdBases
	revDegeneracies = len(reverse) - revBases
	fwdError = fwdDegeneracies / len(forward)
	revError = revDegeneracies / len(forward)

	if (fwdError > revError):
		if (fwdError > errorRate):
			return fwdError
	if (revError > fwdError):
		if (revError > errorRate):
			return revError
	if (fwdError == revError):
		errorRate = fwdError
	return errorRate


# Format all the sequence names in the directory so they're standard: sampleName.readDirection.fastq.gz
# Currently assumes files emerge from the MiSeq as sampleName_Sn_L001_Rn_001.fastq.gz
def renameFiles(directory):
	files = glob.glob(directory + "*.fastq.gz")
	for filename in files:
		path,file = os.path.split(filename)
		filenameSplit = re.split('_',file)
		#print filename,path,file,filenameSplit

		if (len(filenameSplit) == 1):
			return None #break as the files must be formatted already and this is a re-run
		else:
			sampleName = filenameSplit[0]
			readDirection = filenameSplit[3]
			newFilename = path + "/" + sampleName + "." + readDirection + ".fastq.gz"
			rename = os.rename(filename,newFilename)
			print "renamed:\t",filename,"\tto\t",newFilename
	return True

def processArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('--data', help='OPTIONAL: Relative location of the data if not this folder', default=".")
	parser.add_argument('--qiime', help='OPTIONAL: Add this flag if you want to prepare the data for input to QIIME', action='store_true')
	parser.add_argument('--forward', help='Forward primer sequence e.g. ACTGACGTACGT', required=True)
	parser.add_argument('--reverse', help='Reverse primer sequence e.g. ACTGACGTACGT', required=True)
	parser.add_argument('--threads', help='Number of threads to use (bioservc & bioservd max = 80)', default=70, type=int)
	parser.add_argument('--histograms', help='OPTIONAL: Add this flag if you want histograms of the passedQC sequences for each sample', action='store_true')

	args = parser.parse_args()
	return args

if __name__ == '__main__':
	main()
