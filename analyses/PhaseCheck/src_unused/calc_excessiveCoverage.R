# ~~~~ Andrew Frank (c) 2017

setwd("/Users/Andrew/Drive/Work/PlestiodonAnchoredPhylogenomics/PhaseCheck")

tbl <- read.table("./data/samtools_coverage/Combined.coverage")
colnames(tbl) <- c("Locus","Pos","Cov","Sample")

tbl.bysample <- split(tbl, tbl$Sample)
tbl.byboth <- lapply(tbl.bysample, function(x) split(x,x$Locus) )

covDiffs.byboth <- lapply(
  tbl.byboth,
  function(x) {
    lapply(
      x,
      function(y) {
        sample <- y$Sample[1]
        locus <- y$Locus[1]
        pos <- y$Pos[1:(length(y[ ,2])-1)]
        diff <- diff(y$Cov)
        abs_diff <- abs(diff)
        sign <- sign(diff)
        df <- data.frame(sample,locus,pos,abs_diff,sign)
        return(df)
      }
    )
  }
)

covDiffs.combn <- lapply(
  covDiffs.byboth,
  function(x) {
    df <- do.call(rbind,x)
    rownames(df) <- NULL
    return(df)
  }
)
covDiffs.srtd <- lapply(
  covDiffs.combn,
  function(x) {
    pctile.cutoff <- quantile(x$abs_diff, probs = 0.9999)
    cutoff.pos <- which(x$abs_diff >= pctile.cutoff)
    return(x[cutoff.pos, ])
  }
)
