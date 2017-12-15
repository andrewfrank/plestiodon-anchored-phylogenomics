# ~~~~ Andrew Frank (c) 2016

# TO-DO
# Continue working on adding hq.qual column, so far process.samtools_mpileup_vcf
# has been edited, the rest of the functions need to work with it

# 1) Add site distance to end of sequence measure
# 2) Add high quality depth measure (based off DP4)
# 3) Create binomial metric based on depth and allele ratios

source("~/Develop/PhyloFunctions/PhyloFunctions.R")

library(tools)
library(gtools)
library(gdata)

# Functions

process.samtools_mpileup_vcf <- function(file) {

  vcf.table.raw <- try(read.table(  # read in the reduced_representation.vcf
    file = file,                    #   text file
    skip = 29,                      # the first 29 lines contain notes, skip
    colClasses = "character"))      #   those

  # Because some vcf files are empty, they can throw an error - I used the try
  # function above, and if try doesn't capture an error, the function can
  # proceed
  if (! "try-error" %in% class(vcf.table.raw)) {

    # Data clean-up
    locus <- gsub(
      pattern = ".+_L(\\d+).+",
      replacement = "\\1",
      x = file)
    taxon <- gsub(
      pattern = "I\\d+_(.+)_Squamata_.+",
      replacement = "\\1",
      x = vcf.table.raw[ ,1])
    site <- vcf.table.raw[ ,2]

    ref <- vcf.table.raw[ ,4]
    alt <- vcf.table.raw[ ,5]

    qual <- vcf.table.raw[ ,6]
    depth <- gsub(
      pattern = "DP=(\\d+)\\;.+",
      replacement = "\\1",
      x = vcf.table.raw[ ,8])

    allele.depths.raw <- gsub(
      pattern = ".+DP4=(.+)\\;.+",
      replacement = "\\1",
      x = vcf.table.raw[ ,8])
    allele.depths <- strsplit(
      allele.depths.raw,
      ",")
    allele.depths.numeric <- lapply(
      allele.depths,
      as.numeric)
    hq.qual <- lapply(
      allele.depths.numeric,
      sum)

    ref.ratio <- sapply(
      allele.depths,
      function(x) {
        x <- as.numeric(x)
        return(sum(x[1:2])/sum(x[1:4]))
      })
    alt.ratio <- 1 - ref.ratio

    # Write the table
    vcf.table <- data.frame(
      site = as.numeric(site),
      locus = as.numeric(locus),
      taxon = taxon,
      ref = ref,
      alt = alt,
      ref.ratio = ref.ratio,
      alt.ratio = alt.ratio,
      qual = as.numeric(qual),
      hq.qual = hq.qual,
      depth = as.numeric(depth),
      stringsAsFactors = FALSE)

    return(vcf.table)

  } else if ("try-error" %in% class(vcf.table.raw)) {

    return("Empty VCF")             # This is returned if the vcf file is empty

  }

}

process.HapCompass_MWER_solution <- function(file) {

  raw.text <- try(readLines(file))  # read in the MWER_solution.txt file

  length.status <- length(raw.text) != 0
  error.status <- ! "try-error" %in% class(raw.text)

  # Because some MWER-solution files are empty, they can throw an error - I used
  # the try function above, and if try doesn't capture an error, the function
  # can proceed
  if ( all(length.status,error.status) ) {

    # Get table boundaries
    header.lines <- grep(           # find which lines contain the word "BLOCK",
      pattern = "BLOCK",            #   which indicates a table header
      x = raw.text)

    end.lines <- grep(              # find blank lines, which indicate the end
      pattern = "^$",               #   of a table
      x = raw.text)

    table.lines <- data.frame(      # put header lines and end line together for
      header.lines,                 #   use in processing the raw text in the
      end.lines)                    #   function below

    hc.tables.raw <- apply(
      table.lines,
      MARGIN = 1,
      function(x) {
        read.table(
          text = raw.text,
          skip = x[1],
          sep = "\t",
          nrows = x[2] - (x[1] + 1))
      })

    hc.tables.blocked <- lapply(
      seq_along(hc.tables.raw),
      function(x) {
        hc.tables.raw[[x]]$V6 <- x
        return(hc.tables.raw[[x]])
      })

    # Get allele block scores
    split.headers <- sapply(        # read in header lines, split into different
      raw.text[header.lines],       #   components of the header :
      strsplit,                     # <start_position_of_block> <end_position>
      split = "\t")                 #   <start_snp_number> <end_snp_number>
                                    #   <score>

    scores <- sapply(               # get the score from headers
      split.headers,                # a neat one-liner from:
      "[[",                         # http://stackoverflow.com/questions/2803460/how-to-get-the-second-sub-element-of-every-element-in-a-list
      6)

    hc.tables.scored <- lapply(
      seq_along(hc.tables.blocked),
      function(x) {
        hc.tables.blocked[[x]]$V7 <- scores[x]
        return(hc.tables.blocked[[x]])
      })

    hc.table.combn <- do.call(
      rbind,
      hc.tables.scored)

    hc.table <- data.frame(
      site = as.numeric(hc.table.combn[ ,2]),
      block.hc = as.numeric(hc.table.combn[ ,6]),
      block.score.hc = as.numeric(hc.table.combn[ ,7]),
      phase.hc = as.numeric(hc.table.combn[ ,4]),
      stringsAsFactors = FALSE)

    return(hc.table)

  } else if ( ! all(length.status,error.status) ) {

    # This is returned if the MWER file is empty
    return("Empty and/or absent MWER solution")

  }

}

phase.vcf_table <- function(vcf.table,hc.table) {

  if (is.data.frame(hc.table) == TRUE) {

    combn.table <- merge(
      vcf.table,
      hc.table,
      all = TRUE)
    missing.phase <- which(is.na(combn.table$phase.hc))
    combn.table$phase.hc[missing.phase] <- "*"

    block.sizes <- table(combn.table$block.hc)

    if (length(block.sizes) > 1) {

      longest.allele.block <- which(
        block.sizes ==
        max(block.sizes))

      if (length(longest.allele.block) == 1) {

        longest.block.id <- which(
          combn.table$block.hc ==
          longest.allele.block)
        combn.table$phase.hc[-longest.block.id] <- "*"

      } else if (length(longest.allele.block) > 1) {

        block.scores <- unique(combn.table$block.score.hc)

        if (length(block.scores) > 1) {

          highest.score.id <- which(
            combn.table$block.score.hc ==
            max(block.scores, na.rm = TRUE))
          combn.table$phase.hc[-highest.score.id] <- "*"

        } else if (length(block.scores) == 1) {

          qual.avgs <- sapply(
            names(block.sizes),
            function (y) {
              block.ids <- which(
                combn.table$block.hc == y)
              qual.avg <- mean(combn.table$qual[block.ids])})
          highest.qual.block <- which(
            qual.avgs ==
            max(qual.avgs))
          highest.qual.id <- which(
            combn.table$block.hc ==
            highest.qual.block)
          combn.table$phase.hc[-highest.qual.id] <- "*"

        }
      }
    }
  } else if (is.data.frame(hc.table) == FALSE) {

    combn.table <- vcf.table
    combn.table$block.hc <- NA
    combn.table$block.score.hc <- NA
    combn.table$phase.hc <- "*"

  }
  return(combn.table)

}

add.alignment_phase <- function(phased.table,alignment) {

  sites <- phased.table$site
  taxon <- grep(
    pattern = phased.table$taxon[1],
    x = rownames(alignment))

  seq.raw <- alignment[taxon,sites]
  seq <- as.data.frame(
    x = t(toupper(as.character(seq.raw))),
    stringsAsFactors = FALSE)

  phase.ARL <- sapply(
    seq(1,nrow(seq)),
    function (x) {
             if (seq[x,1] == phased.table$ref[x] &&
                 seq[x,2] == phased.table$alt[x]) {
                   phase.ARL <- "0" # matches ref and alt bases
      } else if (seq[x,1] == phased.table$alt[x] &&
                 seq[x,2] == phased.table$ref[x]) {
                   phase.ARL <- "1" # opposite ref and alt bases
      } else if (seq[x,1] == phased.table$ref[x] &&
                 seq[x,2] == phased.table$ref[x]) {
                   phase.ARL <- "R" # homozygous ref base
      } else if (seq[x,1] == phased.table$alt[x] &&
                 seq[x,2] == phased.table$alt[x]) {
                   phase.ARL <- "A" # homozygous alt base
      } else if (seq[x,1] == phased.table$ref[x] &&
                 seq[x,2] != phased.table$alt[x]) {
                   phase.ARL <- "Ra" # ref matches alt different
      } else if (seq[x,1] != phased.table$ref[x] &&
                 seq[x,2] == phased.table$alt[x]) {
                   phase.ARL <- "Ar" # alt matches ref different
      } else if (seq[x,1] != phased.table$ref[x] &&
                 seq[x,2] != phased.table$alt[x]) {
                   phase.ARL <- "*" # ref and alt are ambiguous
      }
    })

  phased.table$phase.ARL <- phase.ARL
  return(phased.table)

}

modifyPhasing.byMismatches <- function(phased.table2) {

  phase.hc.mod <- apply(
    phased.table2,
    MARGIN = 1,
    function(x) {
      if (x[12] == "0" && x[13] == "1") {
        phase.hc.mod <- "*"
      } else if (x[12] == "1" && x[13] == "0") {
        phase.hc.mod <- "*"
      } else phase.hc.mod <- x[12]
    })
  phased.table2$phase.hc.mod <- phase.hc.mod
  return(phased.table2)

}

taxa.data.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/AnchoredPhylo_SamplingScheme.csv"

align.input.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data/Alignments_Raw_phylip_cleaned"
hc.input.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/data/Raw_Reads"

phaseinput.output.path <- "/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/analyses/PHASE/input"

taxa.data.raw <- read.csv(taxa.data.path)

align.files <- list.files(align.input.path, full.names = TRUE)
names(align.files) <- gsub(
  pattern = "T176_L(\\d+).phylip",
  replacement = "\\1",
  x = basename(align.files))
align.files <- align.files[-1]    # remove mtDNA locus

alignments <- lapply(
	align.files,
	read.dna,
	format = "sequential"
)
alignments <- lapply(
  alignments,
  function(x) {
    rownames(x) <- gsub(
      pattern = "\\t",
      replacement = "",
      x = rownames(x))
    return(x)})

# c(1:7,9:13,15,16),
# Here I'll start the loop over each taxon
phased.tables.alltaxa <- lapply(
  c(1:7,9:13,15,16),
  function(i) {

  # Get taxon specific file path for HapCompass folder
  hc.input.path.complete <- file.path(
    hc.input.path,
    taxa.data.raw$LemmonLab_ID[i],
    "anchored_rephasing")

  # Load in the VCF tables - these represent polymorphic sites recognized by
  #   samtools mpileup
  ## A subset of these are empty, either because the alignment never contained
  ##  that taxon in the first place, or no polymorphic sites were recognized
  vcf.files <- list.files(
    path = hc.input.path.complete,
    pattern = "[.vcf]$",
    full.names = TRUE)
  names(vcf.files) <- gsub(
    pattern = "I\\d+.T176_L(\\d+).vcf",
    replacement = "\\1",
    x = basename(vcf.files))

  vcf.tables <- lapply(
    vcf.files,
    process.samtools_mpileup_vcf)

  vcf.tables.cleaned <- vcf.tables[
    which(
      sapply(vcf.tables,class) ==
      "data.frame")]

  # Load in the MWER solution tables - each table represents an allele block that
  #   HapCompass estimated
  ## A subset of these are empty, either because the alignment never contained
  ##  that taxon in the first place, or there were no polymorphic sites recognized
  ##  or because HapCompass couldn't put together allele blocks larger than 1 site
  hc.files <- list.files(
    path = file.path(
      hc.input.path,
      taxa.data.raw$LemmonLab_ID[i],
      "anchored_rephasing"),
    pattern = "_MWER_solution.txt$",
    recursive = TRUE,
    full.names = TRUE)
  names(hc.files) <- gsub(
    pattern = "T176_L(\\d+)_MWER_solution.txt",
    replacement = "\\1",
    x = basename(hc.files))

  hc.tables <- lapply(
    hc.files,
    process.HapCompass_MWER_solution)

  # Make tables with HapCompass phased data, and when there are multiple allele
  # blocks, pick the best allele block (longest, highest score, highest qual)
  phased.tables <- lapply(
    names(vcf.tables.cleaned),
    function (x) {
      phase.vcf_table(
        vcf.table = vcf.tables.cleaned[[x]],
        hc.table = hc.tables[[x]])
      })
  names(phased.tables) <- names(vcf.tables.cleaned)

  # Make tables with ARL phased data
  phased.tables2 <- lapply(
    names(phased.tables),
    function (x) {
      add.alignment_phase(
        phased.table = phased.tables[[x]],
        alignment = alignments[[x]])
    })
  names(phased.tables2) <- names(phased.tables)

  phased.tables2.mod <- lapply(
    phased.tables2,
    modifyPhasing.byMismatches)

  return(phased.tables2.mod)

})

taxon.tables <- lapply(
  phased.tables.alltaxa,
  function (x) do.call(rbind,x))
combn.taxon.table <- do.call(rbind,taxon.tables)
locus.tables <- lapply(
  unique(combn.taxon.table$locus),
  function (x) {
    locus.table <- subset(combn.taxon.table, locus == x)
    locus.table <- locus.table[order(locus.table$site), ]
    rownames(locus.table) <- NULL
    return(locus.table)})
names(locus.tables) <- sapply(
  locus.tables,
  function(x) unique(x$locus))





qual.cutoff <- -10 * log10(0.001)
full.table.cleaned <- full.table.ordered[- which(full.table.ordered$qual < qual.cutoff), ]




test <- sapply(
  locus.tables,
  function(x){

    seq1_nonAmbi.id <- which(
      x$seq1_orig %in%
      c("A","C","G","T"))
    seq2_nonAmbi.id <- which(
      x$seq2_orig %in%
      c("A","C","G","T"))
    nonAmbi.rows <- union(
      seq1_nonAmbi.id,
      seq2_nonAmbi.id)
    ambi.sites.num <- nrow(x) - length(nonAmbi.rows)
    nonAmbi.table <- x[nonAmbi.rows, ]

    seq1_match.id <- which(
      nonAmbi.table$seq1 ==
      nonAmbi.table$seq1_orig)
    seq2_match.id <- which(
      nonAmbi.table$seq2 ==
      nonAmbi.table$seq2_orig)
    matched.rows <- intersect(
      seq1_match.id,
      seq2_match.id)
    matched.sites.num <- length(matched.rows)
    mismatched.table <- nonAmbi.table[-matched.rows, ]

    seq1_seq1orig_match.id <- which(
      mismatched.table$seq1 ==
      mismatched.table$seq1_orig)
    seq1_seq2orig_match.id <- which(
      mismatched.table$seq1 ==
      mismatched.table$seq2_orig)
    seq1_homo.rows <- intersect(
      seq1_seq1orig_match.id,
      seq1_seq2orig_match.id)

    seq2_seq2orig_match.id <- which(
      mismatched.table$seq2 ==
      mismatched.table$seq2_orig)
    seq2_seq1orig_match.id <- which(
      mismatched.table$seq2 ==
      mismatched.table$seq1_orig)
    seq2_homo.rows <- intersect(
      seq2_seq2orig_match.id,
      seq2_seq1orig_match.id)

    homo.rows <- union(seq1_homo.rows,seq2_homo.rows)
    homo.sites.num <- length(homo.rows)
    oddball.table <- mismatched.table[-homo.rows, ]

    unknown.sites.num <- table(x$Known)["*"]

    results <- c(
      unknown.sites.num,
      ambi.sites.num,
      matched.sites.num,
      homo.sites.num,
      nrow(oddball.table))

    return(results)

  })

test.df <- as.data.frame(t(test))
colnames(test.df) <- c(
  "Num_Unknown_Sites",
  "Num_Ambi_Sites",
  "Num_PerfectMatch_Sites",
  "Num_ARLhomozygous_Sites",
  "Num_Oddball_Sites")
table <- test.df[order(test.df$Num_Oddball_Sites, decreasing = TRUE), ]
write.csv(table,"~/Desktop/Loci_Report.csv")


phase.prepped.alignments <- lapply(
  names(locus.tables),
  function(x) {

    align.site.composition <- apply(
      as.character(alignments[[x]]),
      MARGIN = 2,
      table)
    lapply(
      names(align.site.composition),
      function(x) {"-" %in% x})

    allele.variant.sites <- unique(locus.tables[[x]]$base)

    reduced.align <- as.data.frame(
      x = toupper(as.character(
        alignments[[x]][ ,as.numeric(
          )])),
      stringsAsFactors = FALSE)
    colnames(reduced.align) <- unique(locus.tables[[x]]$base)

    reduced.align[] <- lapply(reduced.align, function(x) paste0(x, "0"))

    for (i in 1:nrow(locus.tables[[x]])) {
      rows <- grep(
        pattern = locus.tables[[x]][i,1],
        rownames(reduced.align))
      reduced.align[rows,locus.tables[[x]][i,3]] <- c(
        paste0(locus.tables[[x]][i,4],locus.tables[[x]][i,12]),
        paste0(locus.tables[[x]][i,5],locus.tables[[x]][i,12]))
    }
    reduced.align <- reduced.align[mixedsort(colnames(reduced.align))]
    return(reduced.align)})
names(phase.prepped.alignments) <- names(locus.tables)

lapply(
  names(phase.prepped.alignments),
  function (x) {

    prepped.align <- phase.prepped.alignments[[x]]

    ids <- unique(
      gsub(
        "(.+)_\\w+_seq\\d",
        "\\1",
        rownames(prepped.align)))
    rows <- lapply(
      ids,
      grep,
      x = rownames(prepped.align))
    positions <- colnames(prepped.align)

    prepped.align[] <- lapply(
      prepped.align,
      function(x) {
        gsub(
          pattern = "(.).",
          replacement = "\\1",
          x = x)
        })

    prepped.align[] <- lapply(
      prepped.align,
      function(x) {
        gsub(
          pattern = "-",
          replacement = "?",
          x = x)
        })

    sink(
      file.path(
        phaseinput.output.path,
        paste0("T176_L",x,".inp")))

    cat(length(ids),"\n")
    cat(length(positions),"\n")
    cat("P",positions,"\n")
    cat(rep("S",length(positions)),"\n")

    for (i in 1:length(ids)) {
      cat("#",ids[i],"\n")
      cat(as.character(prepped.align[rows[[i]][1],]),"\n")
      cat(as.character(prepped.align[rows[[i]][2],]),"\n")
    }

    sink()

  })

lapply(
  names(phase.prepped.alignments)[1],
  function (x) {

    prepped.align <- phase.prepped.alignments[[x]]

    prepped.align[] <- lapply(
      prepped.align,
      function(x) {
        gsub(
          pattern = ".(.)",
          replacement = "\\1",
          x = x)
        })

    ids <- unique(
      gsub(
        "(.+)_\\w+_seq\\d",
        "\\1",
        rownames(prepped.align)))
    rows <- lapply(
      ids,
      grep,
      x = rownames(prepped.align))
    kept.rows <- sapply(rows,"[",1)

    sink(
      file.path(
        phaseinput.output.path,
        paste0("T176_L",x,".known.txt")))

    for (i in 1:length(ids)) {
      cat(as.character(prepped.align[kept.rows[[i]],]),"\n")
    }

    sink()

  })


# Pick the longest allele block (like L63)
# Remove allele blocks already perfectly called by Alan (like L241)
# Sets remaining:
#   1) allele blocks containing ambi. not called by Alan, but otherwise match (like L6, L357) - flag these
#   2) allele blocks that disagree on homo vs. hetero calls made by Alan (like L382, L86, L375)
#   3) allele blocks that disagree on which bases belong to which allele (like L234)

# Proposed changes: throw out bases under qual cut-off; if it passes cut-off and
#   Alan called it homozygous, go with heterozygous

# Flag loci with multiple allele blocks

# L45 has weird things happening starting around 475 for all taxa
# How many cases are there when my alternative base disagrees with ARLs?

# EULA_02 L348 had some kind of error, multiple alt bases
