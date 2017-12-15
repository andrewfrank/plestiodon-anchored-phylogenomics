# ~~~~ Andrew Frank (c) 2015

library(ape)
library(gtools)
library(spider)
source("~/Code/Functions.R")

input.path <- "/Users/Andrew/Google Drive/Work/PlestiodonAnchoredPhylogenomics/data/AlignmentsMod_phylip"
output.path <- "/Users/Andrew/Google Drive/Work/PlestiodonAnchoredPhylogenomics/data/AlignmentsMod_phylipReduced"

input.files <- list.files(input.path, full.names = TRUE)
names(input.files) <- list.files(input.path, full.names = FALSE)
input.files <- mixedsort(input.files)

alignments <- lapply(
	input.files,
	read.dna,
	format = "sequential"
)

stripped.alignments <- lapply(
	alignments,
	function(x) {
		missing.data <- checkDNA(x)
		missing.taxa <- which(missing.data == ncol(x))
		stripped.alignment <- x[-missing.taxa,]
		return(stripped.alignment)
	}
)
names(stripped.alignments) <- names(input.files)

lapply(
	names(stripped.alignments),
	function(x) {
		write.dna(
			stripped.alignments[[x]],
			file = file.path(output.path,x),
			format = "sequential",
			nbcol = -1,
			colsep = ""
		)
	}
)