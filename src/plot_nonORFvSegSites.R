# ~~~~ Andrew Frank (c) 2017

library(ape)
library(IRanges)

setwd("~/Drive/Work/PlestiodonAnchoredPhylogenomics")

algn.files <- list.files(
  "./data/Alignments_Raw_phylip_cleaned",
  )
algns <- lapply(
  file.path("./data/Alignments_Raw_phylip_cleaned",algn.files),
  read.dna,
  format = "seq"
  )
names(algns) <- algn.files
algns <- algns[-1]

tbl <- read.table("./PhylogeneticAnalyses/ORFs.bed")
tbl[ ,1] <- gsub("_consensus_sequence","",tbl[ ,1])
tbl.split <- split(tbl,tbl[ ,1])
range.widths <- sapply(
  tbl.split,
  function(x) {
    ir <- IRanges(x[ ,2], x[ ,3])
    ir.r <- reduce(ir)
    w <- ir.r@width
    sum.w <- sum(w)
    return(sum.w)
  }
  )

seq.lengths <- sapply(
  algns,
  ncol
  )
nonORF.sites <- seq.lengths - range.widths
seg.sites <- sapply(
  algns,
  function(x) length(seg.sites(x))
  )
