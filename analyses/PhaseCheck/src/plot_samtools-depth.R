# ~~~~ Andrew Frank (c) 2017

library(ggplot2)
library(ggplus)

setwd("/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/PhaseCheck")

tbl <- read.table("./data/samtools_coverage/Combined.coverage")
colnames(tbl) <- c("Locus","Pos","Cov","Sample")
mean.coverage <- aggregate(tbl$Cov~tbl$Locus, FUN=mean)
tbl$Mean_Cov <- 0

for (i in 1:length(levels(tbl$Locus))) {
  tbl.pos <- which(tbl$Locus == levels(tbl$Locus)[i])          # cells to fill in tbl
  mc.pos <- which(mean.coverage[ ,1] == levels(tbl$Locus)[i])  # row of interest in mean.coverage
  mc.val <- mean.coverage[mc.pos,2]                            # mean coverage value to fill
  tbl$Mean_Cov[tbl.pos] <- mc.val
}

plots <- ggplot(tbl, aes(Pos, Cov)) + geom_line(aes(colour = Locus, linetype = Sample))

pdf("./plots/bwa-mem_coverage.pdf")
plots.pdf <- facet_multiple(plot=plots, facets="Mean_Cov", ncol = 1, nrow = 6)
dev.off()
