# ~~~~ Andrew Frank (c) 2015

source("~/Develop/PhyloFunctions/PhyloFunctions.R")

library(ape)
library(gtools)
library(spider)
library(plyr)

input.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data/Alignments_Raw_phylip_cleaned"
output.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data"

input.files <- list.files(input.path, full.names = TRUE)
names(input.files) <- list.files(input.path, full.names = FALSE)
input.files <- input.files[-1]

alignments <- lapply(
	input.files,
	read.dna,
	format = "sequential"
)

missing.data <- lapply(
	alignments,
	function(x) {
		missing.data <- checkDNA(
			x,
			gapsAsMissing = FALSE)
		names(missing.data) <- gsub(
			pattern = "\\t",
			replacement = "",
			x = names(missing.data))
		percent.missing.data <- round((missing.data / ncol(x)) * 100, digits = 2)
		missing.data.df <- data.frame(
			count.missing.data = missing.data,
			percent.missing.data = percent.missing.data)
		return(missing.data.df)
	})

ambi.data <- lapply(
	alignments,
	function(x) {
		alignment.tables <- apply(
			as.character(x),
			MARGIN = 1,
			table)
		ambi.data <- sapply(
			alignment.tables,
			function (y) {
				ambi.site.id <- which(
					! names(y) %in%
					c("a","c","g","t","n","?","-"))
				return(sum(y[ambi.site.id]))
			})
		names(ambi.data) <- gsub(
			pattern = "\\t",
			replacement = "",
			x = names(ambi.data))
		percent.ambi.data <- round((ambi.data / ncol(x)) * 100, digits = 2)
		ambi.data.df <- data.frame(
			count.ambi.data = ambi.data,
			percent.ambi.data = percent.ambi.data)
		return(ambi.data.df)
	})

locus.data <- lapply(
	names(missing.data),
	function(x) {
		merge(
			missing.data[[x]],
			ambi.data[[x]],
			by = 0)
	})
names(locus.data) <- names(missing.data)

lapply(
	names(locus.data),
	function(x) {
		write.csv(
			file = file.path(output.path,"Alignments_Raw_LocusData",x),
			x = locus.data[[x]])
		})

# Sum percent missing data and percent ambiguous data
locus.data.cat <- rbind.fill(locus.data[-1])
locus.data.totals <- ddply(
	locus.data.cat,
	.(Row.names),
	function(x) colSums(x[,-1], na.rm = TRUE))
locus.data.totals <- locus.data.totals[,c(
	"Row.names",
	"count.missing.data",
	"count.ambi.data")]

write.csv(
	file = file.path(
		output.path,
		"Alignments_Raw_LocusData","T176_AllNuc.csv"),
	x = locus.data.totals)

locus.data.melt <- melt(locus.data.totals)

pdf(
	file = file.path(
		output.path,"MissingData_AmbiData_byTaxa.pdf"),
	width = 11,
	height = 8,
	paper = "USr")
ggplot(
	locus.data.melt,
	aes(x = factor(Row.names), y = value, fill = variable)) +
geom_bar(stat="identity", position = "dodge", color = "black") +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
ggtitle("Total Missing and Ambiguous Sites per Taxa") +
labs(x = "Taxa", y = "Count (bp)") +
scale_fill_manual(
	values = c("black","white"),
	name = "",
	labels = c("?s and Ns", "Ambiguous Bases"))
dev.off()

length.alignments <- sapply(alignments,ncol)

percent.variable.sites <- sapply(
	alignments,
	function(x) {
		x <- del.colgapsonly(x, threshold = 1/nrow(x))
		variable.sites <- seg.sites(x)
		variable.sites.num <- length(variable.sites)
		percent.variable.sites <- variable.sites.num / ncol(x) * 100
	}
)
names(percent.variable.sites) <- names(input.files)

percent.phylo.sites <- sapply(
	alignments,
	function(x) {
		x <- del.colgapsonly(x, threshold = 1/nrow(x))
		phylo.sites <- phylo.sites(x, phased = TRUE)
		phylo.sites.num <- length(phylo.sites)
		percent.phylo.sites <- phylo.sites.num / ncol(x) * 100
	}
)
names(percent.phylo.sites) <- names(input.files)

length.alignments <- sort(length.alignments, decreasing = TRUE)
percent.missing.data.ordered <- sort(percent.missing.data, decreasing = TRUE)

varphylo <- cbind(percent.variable.sites,percent.phylo.sites)
varphylo.df <- as.data.frame(varphylo)
varphylo.df <- arrange(varphylo.df,desc(percent.variable.sites))

barplot(
	length.alignments,
	border = NA,
	xaxt = "n",
	xlab = "Loci",
	ylab = "bp per locus",
	ylim = c(0,2500),
	space = 0
)

barplot(
	percent.missing.data.ordered,
	border = NA,
	xaxt = "n",
	xlab = "Loci",
	ylab = "Percent of missing data per locus",
	ylim = c(0,100),
	space = 0
)

barplot(
	varphylo.df[,1],
	border = NA,
	xaxt = "n",
	xlab = "Loci",
	ylab = "Percent of sites per locus",
	ylim = c(0,30),
	space = 0
)

barplot(
	varphylo.df[,2],
	col = "black",
	border = NA,
	xaxt = "n",
	ylim = c(0,30),
	space = 0,
	add = TRUE
)
