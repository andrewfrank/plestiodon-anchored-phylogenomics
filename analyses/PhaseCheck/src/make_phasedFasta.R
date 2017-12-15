# ~~~~ Andrew Frank (c) 2017

make.ambiguous.base <- function(x) {
  if ( ("A" %in% x) && ("C" %in% x) ) { return("M") }
  if ( ("A" %in% x) && ("G" %in% x) ) { return("R") }
  if ( ("A" %in% x) && ("T" %in% x) ) { return("W") }
  if ( ("C" %in% x) && ("G" %in% x) ) { return("S") }
  if ( ("C" %in% x) && ("T" %in% x) ) { return("Y") }
  if ( ("G" %in% x) && ("T" %in% x) ) { return("K") }
}

library(stringr)

setwd("/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/PhaseCheck")

tbl <- read.table(
  file = "./data/hc_tables/I7487_srtdhc.table",
  header = TRUE,
  na.strings = "NA",
  stringsAsFactors = FALSE)
split.gt <- str_split_fixed(tbl[ ,5], "/", 2)
tbl$SEQA <- split.gt[ ,1]
tbl$SEQB <- split.gt[ ,2]

tbl.inv <- tbl[which(tbl$SEQA == tbl$SEQB), ]
tbl.var <- tbl[which(tbl$SEQA != tbl$SEQB), ]
tbl.var.uph <- tbl.var[which(is.na(tbl.var[ ,6])), ]

SEQA.len <- unlist(lapply(tbl$SEQA,nchar))
SEQB.len <- unlist(lapply(tbl$SEQB,nchar))
indels <- union(
  which(SEQA.len != 1),
  which(SEQB.len != 1))

tbl.indels <-tbl[indels, ]
tbl.snps


unphased.SEQA.len <- unlist(lapply(tbl.unphased$SEQA,nchar))
unphased.SEQB.len <- unlist(lapply(tbl.unphased$SEQB,nchar))
unphased.INDELs <- union(
  which(unphased.SEQA.len != 1),
  which(unphased.SEQB.len != 1))
tbl.unphased.INDELs <- tbl.unphased[unphased.INDELs, ]
tbl.unphased.SNPs <- tbl.unphased[-unphased.INDELs, ]

ambiguous.SNPs <- unlist(
  apply(
    tbl.unphased.SNPs,
    MARGIN = 1,
    function(x) make.ambiguous.base(c(x[7], x[8]))
  )
)
tbl.unphased.SNPs$SEQA <- ambiguous.SNPs
tbl.unphased.SNPs$SEQB <- ambiguous.SNPs
