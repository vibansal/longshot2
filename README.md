
# Introduction

This repository contains Longshot2 (implemented in python) that implements a two-pass approach to variant calling (both SNVs and short indels) for long reads. Longshot2 uses a C++ implementation of the POA algorithm (https://github.com/rvaser/spoa) and the Longshot tool (github.com/pjedge/longshot). 

Longshot is a variant calling tool for diploid genomes using long error prone reads such as Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT)
available from https://github.com/pjedge/longshot

## Algorithm overview 

1. In the first pass, SNVs are called using Longshot and reads are separated by haplotype 
2. POA-based consensus algorithm is applied to the haplotype-separated reads to identify variants (both SNVs and indels)
3. POA-based variants are merged with the SNVs from step (1) using bcftools 
4. The merged set of variants are used as input for a second run of Longshot

# Installation 

1. Download and compile the python wrapper for C++ spoa library:

```
git clone https://github.com/vibansal/spoa_python.git
cd spoa_python
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
python3.x setup.py build_ext --inplace
cd ..
cp -r spoa_python/spoa .
```

2. Install longshot (https://github.com/pjedge/longshot)

3. Install other dependencies: bgzip, tabix, bcftools 

4. edit the config.txt to change the path to longshot and bcftools

# Usage 

python3.6 longshot2.py [args] 

The python program expects the same arguments as used for longshot. Currently, the "--region" option is required.

Example command:

```
python3.6 longshot2.py --bam HG002.longreads.bam --region chr20 --ref GRCh38.fasta  --out HG002.chr20.longshot2.vcf 
```


