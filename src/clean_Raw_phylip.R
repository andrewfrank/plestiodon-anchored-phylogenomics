# ~~~~ Andrew Frank (c) 2015

source("~/Develop/PhyloFunctions/PhyloFunctions.R")

library(ape)

input.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data/Alignments_Raw_phylip"
output.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data/Alignments_Raw_phylip_cleaned"

input.files <- list.files(
	input.path,
	pattern = "^T176_L\\d+.phylip$",
	full.names = TRUE)
names(input.files) <- basename(input.files)

alignments <- lapply(
	input.files,
	read.dna,
	format = "sequential")

alignments.renamed <- lapply(
  alignments,
  clean.names)

lapply(
  names(alignments.renamed),
  function(x) {
    write.dna(
      file = file.path(output.path,x),
      x = alignments.renamed[[x]],
      format = "sequential",
      nbcol = -1,
      colsep = "")
  })
