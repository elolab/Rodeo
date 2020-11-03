# Rodeo description

Rodeo (RObust DEcOnvolution) is a robust deconvolution tool that provides estimates of how strongly each cell type expresses each gene based on bulk expression data and cell type proportions. The current version 1.0 of the package includes one function intended for external usage: Rodeo. Please contact us (maria.jaakkola@utu.fi) in case of bugs, missing documentation, or problems with installation. 

### Installation and usage

Besides R package Rodeo (available also at https://elolab.utu.fi/software/), also package MASS should be installed (comes with R by default) and activated with R command library(MASS) before using Rodeo. Rodeo can be installed by opening R and typing devtools::install_github("elolab/Rodeo") (requires package devtools to be installed).


Below is a code example of how to run Rodeo after it has been installed. The example is very simplistic and can not be directly copy pasted if e.g. the input data is stored differently.

	`\# Load required packages`
	`library(MASS)`
	`library(Rodeo)`
	
	`\# Read Input data`
	`setwd("path/to/my/input/files")`
	`E = read.table("BulkData.txt", header=T, row.names=1)`
	`C = read.table("CellTypeProportions.txt", header=T, row.names=1)`
	
	`\# Run Rodeo and save the results`
	`S = Rodeo(E, C)`
	`write.table(S, file="EstimatedS.txt", sep="\\t", quote=F)`


### Input and output

The current version of Rodeo has two input parameters: 
| Input      | Description |
| ----------- | ----------- |
| E      | a bulk expression matrix with named rows (genes) and columns (samples)       |
| C   | a cell type proportion matrix with named rows (cell types) and columns (samples) for the samples in **E**. Each sample (i.e. column) should sum to 1.        |


The main function (and the only one meant for external use) Rodeo returns a matrix S where columns are cell types, rows are genes, and elements describe how strongly the cell type expresses the gene. 