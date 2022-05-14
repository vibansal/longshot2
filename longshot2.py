#!/usr/bin/env python3
import os
import time
import argparse
import sys
import subprocess
import poacaller
import pysam
import gzip
import multiprocessing

def read_config(args,filename='config.txt'):
	#args = {}
	File = open(filename)
	for line in File:
		option = line.strip().split('=')
		args[option[0]] = option[1]
	File.close()
	#return args

def extract_args_longshot(args_list):
	args_table= {}
	include_vec = [1]*len(args_list)
	for i in range(len(args_list)-1):
		if args_list[i] == '--region'  or args_list[i] == '-r': 
			args_table['region'] = args_list[i+1]
			include_vec[i] = 0
			include_vec[i+1]=0
		elif args_list[i] == '--ref'  or args_list[i] == '-f': args_table['fasta'] = args_list[i+1]
		elif args_list[i] == '--out'  or args_list[i] == '-o': 
			args_table['out'] = args_list[i+1]
			include_vec[i] = 0
			include_vec[i+1]=0
		elif args_list[i] == '--out_bam'  or args_list[i] == '-O': 
			args_table['haplobam'] = args_list[i+1]
			include_vec[i] = 0
			include_vec[i+1]=0
		elif args_list[i] == '--bam'  or args_list[i] == '-b': args_table['bam'] = args_list[i+1]
		elif args_list[i] == '-F': args_table['F'] = True
	filtered_args = [args_list[i] for i in range(len(args_list)) if include_vec[i] ==1 ]
	#print(filtered_args,'filter'); sys.exit()
	return args_table,filtered_args


def call_POA(args_table,bam,fasta,outdir,inputvcf,region):
	parser = poacaller.parseargs(empty=True)
	poa_args=parser.parse_args([]) ## empty 
	poa_args.bam = bam
	poa_args.ref = fasta
	poa_args.out = outdir + '/poa'
	poa_args.vcf = inputvcf ## longshot snvs vcf file
	poa_args.region = region
	#print(poa_args)

	print('\n\nrunning POA variant detection on haplotagged bam file',poa_args.bam,'in region',region,file=sys.stderr)
	varcaller = poacaller.VarCall(poa_args)
	varcaller.parse_region(region)
	varcaller.process_bam()
	varcaller.close_filehandles()
	poacaller.combinevcfs(poa_args.out,poa_args.ref,BCFTOOLS=args_table['BCFTOOLS'])

	return poa_args.out +  '.final.diploid.vcf.gz' ## output vcf file


#bcftools filter -e 'FILTER~"dn" & INFO/SOURCE = "INP"' poa-test/poa.final.diploid.vcf.gz 
def filter_dense_variants(filename,outfile):
	of = open(outfile,'w')
	of_lg = open(outfile+'.largeindels','w')
	if filename.endswith(".gz"): f = gzip.open(filename)
	else: f = open(filename)
	for line in f:
		try:  sline=line.decode()
		except AttributeError: sline = line
		if sline[0] == '#': 
			print(sline.strip(),file=of)
			print(sline.strip(),file=of_lg)
		else: 
			var = sline.strip().split('\t')
			indel_length = max(len(var[3]),len(var[4]))
			if indel_length > 50: 
				print(sline.strip(),file=of_lg)
				continue ## filter large indels
			if ('dn' in var[6] and var[7] == 'SOURCE=INP') or 'N' in var[3]: 
				#print('filtering variant',var,file=sys.stdout)
				continue
			elif var[7] == 'SOURCE=INP': 
				print('\t'.join(var[0:10]+['0/1:.']),file=of)
			elif 'dn' not in var[6]: 
				print('\t'.join(var),file=of)
			else: 
				#print('changing INFO field to PASS',var,file=sys.stdout)
				print('\t'.join(var[0:6]+['PASS'] + var[7:]),file=of) 
	f.close()
	of.close()
	of_lg.close()
	

def main():
	args_table,filtered_args = extract_args_longshot(sys.argv[1:])
	read_config(args_table)
	#if 'F' in args_table and os.path.isdir(args_table['out']): 
	if 'bam' not in args_table or 'fasta' not in args_table or os.path.isdir(args_table['OUTDIR']): 
		print('check input arguments',file=sys.stderr)
		sys.exit()
	print(args_table,'\n')

	os.makedirs(args_table['OUTDIR'])

	regions = []
	if 'region' not in args_table: 
		pyfasta = pysam.Fastafile(args_table['fasta'])
		for contig in pyfasta.references: regions.append(contig)
	else: regions.append(args_table['region'])
	
	for region in regions:
		print('\n\nrunning longshot first pass for SNV calling',region,file=sys.stderr)
		args_table['haplobam'] = args_table['OUTDIR'] + '/haplotag.bam'
		args_table['poavcf'] = args_table['OUTDIR'] + '/poa_filtered.vcf'
		args_table['poavcf_gz'] = args_table['OUTDIR'] + '/poa_filtered.vcf.gz'
		args_table['outvcf1'] = args_table['OUTDIR'] + '/snvs.vcf'
		args_table['outvcf2'] = args_table['OUTDIR'] + '/snvs_indels.vcf'
		args_table['region'] = region

		longshot_cmd =args_table['LONGSHOT'] + ' ' + ' '.join(filtered_args) + ' --out_bam ' + args_table['haplobam'] + ' --out ' + args_table['outvcf1'] + ' --region ' + region
		#print(longshot_cmd)
		subprocess.call(longshot_cmd,shell=True)
		subprocess.call('bgzip -f ' + args_table['outvcf1'],shell=True)
		subprocess.call(args_table['TABIX'] + ' -f ' + args_table['outvcf1'] + '.gz',shell=True)

		if not os.path.isfile(args_table['haplobam']): 
			print('no haplotype-tag bam file',file=sys.stderr)
			continue

		output_vcf_poa = call_POA(args_table,args_table['haplobam'],args_table['fasta'],args_table['OUTDIR'],args_table['outvcf1'] + '.gz',region)

		filter_dense_variants(output_vcf_poa,args_table['poavcf'])
		subprocess.call('bgzip -f ' + args_table['poavcf'],shell=True)
		subprocess.call(args_table['TABIX'] + ' -f ' + args_table['poavcf_gz'],shell=True)

		print('\n\nre-running longshot with POA variants as input',file=sys.stderr)
		longshot_cmd =args_table['LONGSHOT'] + ' ' + ' '.join(filtered_args) + ' -v ' + args_table['poavcf_gz']  + ' --out ' + args_table['outvcf2'] + ' --region ' + region + ' -D 100:500:50'
		subprocess.call(longshot_cmd,shell=True)

	subprocess.call('cp ' + args_table['outvcf2'] + ' ' + args_table['out'],shell=True)


if __name__ == "__main__":
	main()
