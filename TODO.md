# build script to take Alan’s alignments and transform them into vcf files

this only detects SNPs from multiple sequence alignments in fasta format (I need indels as well) https://github.com/sanger-pathogens/snp-sites

then use vcftools vcf-merge: http://vcftools.sourceforge.net/perl_module.html#vcf-merge

think about options to take multiple singletons resulting from INDELS (i.e. A,A and C,- and G,-) and merging them together into a single INDEL (A,ACG)

# compare my phased variants to Alan’s phase variants

will need to compare vcf files, see which variants fall under either "shared", "jockusch-only" or "lemmon-only"

# build script which takes phased variants and transforms them into completed alignments

try https://github.com/Bahler-Lab/alignment-from-vcf
but will need to capture indels as well
