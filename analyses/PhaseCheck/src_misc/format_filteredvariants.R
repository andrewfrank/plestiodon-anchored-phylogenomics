# ~~~~ Andrew Frank (c) 2017

setwd("~/Drive/Work/Projects/plestiodon-anchored-phylogenomics/analyses/PhaseCheck")

dat <- read.table(
  file = "./data/variants/VariantCalls_fltrd.vcf",
  skip = 434
)
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","I18066","I18068","I18069","I7487","I7488","I7489","I7490","I7491","I7492","I7493","I7495","I7496","I7497","I7498")

dat.fltrd <- dat[which(dat$FILTER != "PASS"), ]
dat.fltrd.haps <- apply(
  X = dat.fltrd[ ,10:23],
  MARGIN = 2,
  FUN = gsub,
  pattern = "(\\S\\/\\S)\\S+",
  replacement = "\\1"
)
dat.fltrd.hapspresent <- apply(
  X = dat.fltrd.haps,
  MARGIN = 1,
  FUN = function (x) 14 - length(which(x == "0/0")) - length(which(x == "./."))
)
dat.new <- data.frame(
  dat.fltrd[ ,1:9],
  dat.fltrd.haps,
  PRESENCE = dat.fltrd.hapspresent
)

write.table(
  dat.new,
  file = "./docs/VariantCalls_fltrd_modified.tsv",
  quote=FALSE,
  sep='\t',
  col.names = NA
)
