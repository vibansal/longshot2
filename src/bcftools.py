from __future__ import print_function
import os,sys,glob,math,argparse
import subprocess


## -a split-MNPs useful for merging with longshot variants...
def combinevcfs(inputfilename,fasta,BCFTOOLS,INPUT=True,MAX_LENGTH=100,NORMALIZE=True):
	array = ['POA_0','POA_1','POA_2']
	if INPUT: array.append('INP')

	for hap in array:
		outfile=inputfilename + '.' + hap + '.vcf.gz'
		command= BCFTOOLS + ' filter -i \'SOURCE=\"' + hap  + '\"\' ' + inputfilename + ' | ' + BCFTOOLS + ' norm -a ' + ' -f ' + fasta + ' -d exact -O z -o ' + outfile 
		print(command,file=sys.stderr)
		subprocess.call(command,shell=True)
		subprocess.call(BCFTOOLS + ' index ' + outfile,shell=True)

	files = ' '.join([inputfilename + '.' + val + '.vcf.gz' for val in array])
	command=BCFTOOLS + ' merge -m all ' + files + ' --force-samples -i SOURCE:join -o ' + inputfilename + '.final'
	print(command,file=sys.stderr)
	subprocess.call(command,shell=True)
	for hap in array: subprocess.call('rm -rf ' + inputfilename + '.' + hap + '.vcf*',shell=True)
	merge_genotypes(inputfilename + '.final')


## calculate distributions of ratios for snp/ins/del and filter based on that...

def merge_genotypes(input_vcf_file,out_vcf_file=None,SAMPLE_ID='SAMPLE'):  ## multiple genotypes (POA)
	if out_vcf_file == None: out_vcf_file = input_vcf_file + '.diploid.vcf'
	outfile = open(out_vcf_file,'w')
	File = open(input_vcf_file)
	for line in File:
		if line[0] == '#' and line[1] == 'C': 
			variant = line.strip().split('\t')
			for i in range(9): print(variant[i],file=outfile,end='\t')
			print(SAMPLE_ID,file=outfile)
		elif line[0] == '#': 
			print(line.strip(),file=outfile)
		else:
			variant = line.strip().split('\t')
			alleles = len(variant[4].split(','))
			g0 = variant[9].split(':')[0] ## diploid
			g1 = variant[10].split(':')[0]
			g2 = variant[11].split(':')[0]
			if g1 == './.' and g2 == './.': ## only joint POA call
				#frac = float(variant[9].split(':')[5])
				#if frac > 0.75: genotype= '1/1' + variant[9][3:]
				#else: genotype = '0/1' + variant[9][3:]
				genotype = variant[9]
				## 0/1:0:105190-105490:26:26:0.87,0.003,1,1D
			elif g1 != './.' and g2 != './.': ## both haplotypes
				if  g1.split('|')[0] != '0': allele1 = g1.split('|')[0]
				elif  g1.split('|')[1] != '0': allele1 = g1.split('|')[1]
				if  g2.split('|')[0] != '0': allele2 = g2.split('|')[0]
				elif  g2.split('|')[1] != '0': allele2 = g2.split('|')[1]
				genotype = allele1 + '/' + allele2 + variant[10][3:]
			elif g1 == './.' and g2 != './.' and alleles ==1: genotype = variant[11]
			elif g2 == './.' and g1 != './.' and alleles ==1: genotype = variant[10]
			elif g1 == './.' and g2 != './.' and alleles > 1: 
				genotype = '1/2' + variant[11][3:]
			elif g2 == './.' and g1 != './.' and alleles > 1: 
				genotype = '1/2' + variant[10][3:]

			for i in range(9): print(variant[i],file=outfile,end='\t')
			print(genotype,file=outfile)
				
	File.close()
	outfile.close()
	subprocess.call('bgzip ' +  out_vcf_file,shell=True)
	subprocess.call('tabix  ' +  out_vcf_file+'.gz',shell=True)

if __name__ == '__main__':
	## vcf file to process and fasta file
	combinevcfs(sys.argv[1],sys.argv[2])
	
	
