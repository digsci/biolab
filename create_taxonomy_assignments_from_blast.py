import os,sys,re,argparse
from collections import defaultdict

def main():
	options = processArguments()
	taxonomyDB = createTaxonomyDB(options)
	assignedTaxonomy = assignTaxonomy(taxonomyDB,options)

def assignTaxonomy(taxonomyDB,options):
	assigned = {}
	outputfile = open(options.output,"w")
	for line in open(options.blast,'rU'):
		if "pid" not in line:
			linelist = re.split('\t',line)
			otu = linelist[0]
			percid = float(linelist[2])
			evalue = linelist[10]
			hit = linelist[1]
			if (percid >= options.percid):
				taxonomy = taxonomyDB[hit]
				outputline = str(otu) + "\t" + str(taxonomy) + "\t" + str(evalue) + "\t" + str(hit) + "\n"
				outputfile.write(outputline)
			else:
				outputline = str(otu) + "\tNo_blast_hit;\tNone\tNone\n"
				outputfile.write(outputline)
	return "done"

def createTaxonomyDB(options):
	taxonomy = {}
	for line in open(options.taxonomy,'rU'):
		linelist = re.split('\t',line)
		taxonomy[linelist[0]] = linelist[1].rstrip()
	return taxonomy

def processArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('--taxonomy', help='Diatom database taxonomy file (diatoms.taxonomy.FINAL.txt)', required=True)
	parser.add_argument('--percid', help='Percentage ID to filter',required=True, type=float)
	parser.add_argument('--blast', help='Blast results', required=True)
	parser.add_argument('--output', help="Name of the output file", required=True)
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	main()
