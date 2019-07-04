#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
	stop("Usage: Rscript darleq3process.R inputfile outputfile\nExample: Rscript darleq3process.R darleq3-input.xlsx darleq3-output.xlsx")
} else {
	library("darleq3")
	inputfile <- args[1]
	outputfile <- args[2]
	dataframe <- read_DARLEQ(inputfile)
	results <- calc_Metric_EQR(dataframe, metrics="TDI5NGS")
	save_DARLEQ(results, outFile=outputfile)
}

