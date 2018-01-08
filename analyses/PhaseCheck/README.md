# Plestiodon Anchored Phylogenomics: PhaseCheck

**Pipeline order:**
```
run_concatenate.sh
run_sickle.sh
run_mapReads.sh
      |       \
      |        \
      |   run_samtools-depth.sh
      |   plot_samtools-depth.R
      |   build_BEDfilter
      |        /
      |       /       
run_readFilter.sh
run_gatk-HaplotypeCaller.sh
run_finalizeSNPs.sh
run_vcfToPhasedFasta.sh
```

## Read preparation

### Read renaming

I manually reduced the original file names to more easily readable names. Example:

| Original File Name                        | New File Name      |
| P0077_EJ_I7487_GGTAGTTG_L002_R1_001.fastq | I7487_R1_001.fastq |
| P0077_EJ_I7487_GGTAGTTG_L002_R1_002.fastq | I7487_R1_002.fastq |
| P0077_EJ_I7487_GGTAGTTG_L002_R2_001.fastq | I7487_R2_001.fastq |
| P0077_EJ_I7487_GGTAGTTG_L002_R2_002.fastq | I7487_R2_002.fastq |

These are saved in /data/original_reads

### Read file concatenation

run_concatenate.sh simply takes the "001" and "002" read files and combines them into a single fastq.gz file. I.e. I7487_R1_cast.fastq.gz.

Concatenated read files are saved in /data/concat_reads

### Trimming with sickle

run_sickle.sh runs the trimming program sickle to remove adapters and low quality sequence.

Trimmed read files are saved in /data/trimmed_reads

## Read mapping with BWA-MEM

In order to test the phasing done by Alan's program, I need to map the raw reads back to his sequences. I decided that this was the best approach, given that denovo assembly wouldn't create a 1-1 comparison with Alan's sequences. bwa-mem offered a flexible solution between end-to-end and local alignment appropriate for very short reference sequences (confirmed by Jill, so contact her if need help justifying this).

### Index construction

To get sequences to map back to, I created a single strict consensus sequence per locus alignment given to me by Alan. I forced the insertion of a the most common base at a given site in the consensus sequence instead of leaving it ambiguous. This is because bwa-mem randomly assigns bases whenever it sees an N. I'm not sure how it deals with other ambiguous characters.) I used Geneious to output the consensus sequence for each alignment into a fasta file.

The original fasta files for the phased alignments I received from Alan are saved in PlestiodonAnchoredPhylogenomics/data/original_seqs (the overall project directory)

The consensus sequences per alignment are saved in PlestiodonAnchoredPhylogenomics/data/original_seqs (the overall project directory)

I then concatenated these fasta files into a single fasta file, ConsensusSeqs_Raw.fasta:
```
cat *.fasta > ConsensusSeqs_Raw.fasta
```

I saved ConsensusSeqs_Raw.fasta in /data/bwa-index.

There were a limited number of cases where Geneious was unable to call a non-ambiguous base at a given site, because the site was evenly split between two alleles (i.e. 7 Gs and 7 Ts). This resulted in insertion of ambiguous bases into these sites. I found out during the downstream haplotype construction phase (post variant calling and phasing) that having these ambiguous bases persist resulted in bwa-mem and all downstream programs to ignore these sites.

To resolve these ambiguous sites, I did a simple replacement of all ambiguous sites to one of the standard sites. This guarantees bwa-mem will see the site, map accordingly, and then variation at these sites (which I expect), will be detected by HaplotypeCaller and GenotypeGVCFs.

```
sed \
  -e s/M/A/g \
  -e s/R/A/g \
  -e s/W/A/g \
  -e s/S/C/g \
  -e s/Y/C/g \
  -e s/K/G/g \
  -e s/V/A/g \
  -e s/H/A/g \
  -e s/D/A/g \
  -e s/B/C/g \
  concat-consensus_seqs.fasta > concat-consensus_seqs_noAmbi.fasta
```

I then used BWA-MEM to construct an index:
```
bwa index -p ConsensusSeqs -a is ConsensusSeqs.fasta
```

These bwa index files were also saved in /data/bwa-index

***Note***: I created an alternative approach to reference sequence mapping, saved as run_bwa-index.sh in /src_unused. This approach created an index per sample, where each sequence was that sample's particular sequence in an alignment. This lead to downstream issues, since not all samples were present in every alignment from Alan. It was easier to assess each sample across every loci using the entire alignment consensus. These per sample indices are saved in /data_unused/bwa-index

### Create fasta dict with picard's CreateSequenceDictionary

From https://software.broadinstitute.org/gatk/documentation/article?id=1601
```
java -jar $PICARD \
  CreateSequenceDictionary \
  R=../consensus_seqs/concat-consensus_seqs_noAmbi.fasta \
  O=../consensus_seqs/concat-consensus_seqs_noAmbi.dict
```

### Create fasta index with samtools

From https://software.broadinstitute.org/gatk/documentation/article?id=1601

```
samtools faidx ConsensusSeqs.fasta
```

### Read mapping

run_bwa-mem runs bwa-mem on each species' paired end read files, then uses picard to make a BAM index file. The resulting BAM and BAI files are saved in /data/bwa-mappings

### Add read group identifiers with picard's AddOrReplaceReadGroups

### Mark PCR and optical duplicates with picard's MarkDuplicates

## Get per base coverage with samtools depth

To assess bwa-mem's mapping success, I wanted to visualize per-base coverage across each locus for each sample. To obtain per-base coverage, run_samtools-depth gets coverage data. These coverage text files are saved in /data/samtools_coverage

## Search for loci with excessive coverage

To examine for the presence of excess coverage, I visualized coverage data into a series of plots per loci in plot_samtools-depth.R. I saved resulting plots, sorted from lowest average coverage to highest average coverage, in /plots/bwa-mem_coverage.pdf

I noted loci that had plots with visually striking excess coverage. Notably, these loci had some of the highest average coverage. I then uploaded the bwa-mem bam files into Geneious, and examined the read pileups of these noted loci. These loci either contained repeat regions with excessive coverage, or flanking regions with excessive coverage. Many reads that overlapped these regions were soft-clipped by bwa-mem outside the excessive coverage region. This indicates that bwa-mem matched only part of the read to the excessive region, but the rest of the read did not match. This is evidence that these reads are incorrectly mapped. Notably, bwa-mem's soft clipping indicated clear ranges for the excessive coverage regions. I annotated these regions on the sequences in Geneious, and wrote a separate file recording excess coverage regions: /data/excessive_coverage_loci.csv

I then manually constructed a BED file which indicates areas of non-excessive coverage per locus. For most of the locus, this is all of its bases (0 - end of sequence), but for loci with excessive coverage, it has the excessive coverage regions removed (this is essentially the opposite file of excessive_coverage_loci.csv).

***Note***: I attempted to *a priori* assess areas of excess coverage with the script calc_excessiveCoverage.R, saved in /src_unused. The script attempts to detect when, in a given vector of per-base coverage differences across a locus, the greatest changes in coverage occur. However, I had substantial difficulty creating a cutoff for degree in coverage change that accurately reflected what I observed in the coverage plots and the pileups visualized in Geneious.

## Checking called variants against Alan's called variants

```
module load gatk

java -jar $GATK \
  -T ValidateVariants \
  -R ../consensus_seqs/concat-consensus_seqs_noAmbi.fasta \
  -V I7487_lemmons.vcf

java -jar $GATK \
  -T LeftAlignAndTrimVariants \
  -R ../consensus_seqs/concat-consensus_seqs_noAmbi.fasta \
  --variant I7487_lemmons.vcf \
  -o test.vcf

```

```
vcf-merge I7487_lemmons.vcf.gz I7488_lemmons.vcf.gz I7489_lemmons.vcf.gz I7490_lemmons.vcf.gz I7491_lemmons.vcf.gz I7492_lemmons.vcf.gz I7493_lemmons.vcf.gz I7495_lemmons.vcf.gz I7496_lemmons.vcf.gz I7497_lemmons.vcf.gz I7498_lemmons.vcf.gz I18066_lemmons.vcf.gz I18068_lemmons.vcf.gz I18069_lemmons.vcf.gz -c any -t

java -jar $GATK \
  -T CombineVariants \
  -R ../consensus_seqs/concat-consensus_seqs_noAmbi.fasta \
  --variant I7487_lemmons.vcf \
  --variant I7488_lemmons.vcf \
  --variant I7489_lemmons.vcf \
  --variant I7490_lemmons.vcf \
  --variant I7491_lemmons.vcf \
  --variant I7492_lemmons.vcf \
  --variant I7493_lemmons.vcf \
  --variant I7495_lemmons.vcf \
  --variant I7496_lemmons.vcf \
  --variant I7497_lemmons.vcf \
  --variant I7498_lemmons.vcf \
  --variant I18066_lemmons.vcf \
  --variant I18068_lemmons.vcf \
  --variant I18069_lemmons.vcf \
  -genotypeMergeOptions UNIQUIFY \
  -o
  
```
## Filtering reads

I wanted to apply 2 filter to the bwa-mem mapped reads:
1) Remove reads that were internally soft-clipped, which was indicative of off target reads that had erroneously mapped, while keeping reads that were soft-clipped on the edges of the reference sequence, in order to preserve coverage on the edges of the original sequences.
2) Remove read pairs that were exclusively mapped to regions of excess coverage, since they lacked evidence of actually being attached to that particular reference sequence.

I made run_readFilter.sh to perform these two filter steps. It first filters for reads only in primary alignments and are properly mapped (both the read pairs are mapped to the same reference seq). Then it uses awk to filter out all internally soft-clipped reads, leaving only soft-clipped reads at the ends of the ref sequence. It then uses samtools to filter mapped reads by a BED file I generated (saved in /data/excessive_coverage_filter.bed), which saves all reads overlapping each loci's non-excess coverage region. Reads filtered out at this step are saved to a separate file, and tossed reads that have mates in the non-excess coverage region are reincorporated via join. I perform a final samtools sort and index the BAM file with picard.

## Recalibrating bases with GATK's BaseRecalibator

Optional step, not sure if I'll do it.

https://software.broadinstitute.org/gatk/documentation/article?id=1247

## Call variants per-individuals with GATK's HaplotypeCaller

Because

## Call variants / phase across individuals, and filter SNPs & INDELs

### GATK's GenotypeGVCFs

### Filter SNPs and INDELs through GATK's recommended filters

## Split VCF by individual, and output haplotype sequences
