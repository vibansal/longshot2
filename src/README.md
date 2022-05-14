## POA based variant identification

Python code to identify candidate variants from haplotype-separated long read alignments (BAM/CRAM files) for diploid genomes.


## requirements

* pysam
* bgzip 
* tabix 
* bcftools >= v1.12
* parasail


## usage 

```
usage: main.py [-h] [-w [WINDOWSIZE]] [-m [MINREADS]] [--minMQ [MINMQ]]
               [--maxdepth [MAXDEPTH]] [--flank [FLANK]] [--ref [REF]]
               [-b [BAM]] [--no-hap] [--addref] [--usecigar] [--cluster]
               [--compare] [--filterreads] [--pre] [-o [OUT]] [--vcf [VCF]]
               [-i [REGION]] [--bed [BED]]

optional arguments:
  -h, --help            show this help message and exit
  -w [WINDOWSIZE], --windowsize [WINDOWSIZE]
                        window size for POA calculation
  -m [MINREADS], --minreads [MINREADS]
                        minumum reads required for POA consensus
  --minMQ [MINMQ]       minumum mapping quality of reads used for POA
                        consensus
  --maxdepth [MAXDEPTH]
                        maximum number of reads allowed for POA consensus
  --flank [FLANK]       flanking sequence length for consensus, default 30 bp
  --ref [REF], --fasta [REF]
                        reference fasta file to which reads are mapped
  -b [BAM], --bam [BAM]
                        bam file
  --no-hap              do not use haplotype tags from BAM file (if available)
                        for POA
  --addref              add reference sequence as first sequence in POA
  --usecigar            add first read using cigar alignment in POA
  --cluster             group identical reads to speed up processing, for CCS
                        reads
  --compare             compare read counts supporting variant between two
                        haplotypes
  --filterreads         filter out sub-reads that match reference sequence
  --pre                 calculate per-read-error rates for POA ordering
  -o [OUT], --out [OUT]
                        output vcf file,combined
  --vcf [VCF]           VCF file with candidate variants, gzipped and indexed
  -i [REGION], --region [REGION]
                        target genomic interval
  --bed [BED]           bedfile with intervals



```
