
# Introduction

Longshot is a variant calling tool for diploid genomes using long error prone reads such as Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT)
available from https://github.com/pjedge/longshot

This repository contains Longshot2 (implemented in python) that implements a two-pass approach to variant calling (both SNVs and short indels) for long reads. Longshot2 uses the C++ implementation of the POA algorithm (https://github.com/rvaser/spoa) to generate consensus sequences and identify variants. 


# Installation 

1. The first step is to compile the python wrapper for C++ spoa library:

```
cd spoa_python
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
python3 setup.py build_ext --inplace
cd ..
cp -r spoa_python/spoa .
```

2. Install longshot (https://github.com/pjedge/longshot)

3. Install other dependencies: bgzip, tabix, bcftools 

4. edit the config.txt to change the path to longshot and bcftools

# Usage 

python3.6 longshot2.py [args] 

Provide the same arguments as used for longshot. 

