#!/usr/bin/env python3
import os
import time
import argparse
import sys
import subprocess
import poacaller
import pysam
import gzip
from multithreading import *
from concurrent.futures import ThreadPoolExecutor

NCORES=6

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


def call_POA_direct(region,args_table):
	parser = poacaller.parseargs(empty=True)
	poa_args=parser.parse_args([]) ## empty 
	poa_args.bam = args_table[(region,'haplobam')]
	if not os.path.isfile(poa_args.bam): 
		print('no haplotype-tag bam file',file=sys.stderr)
		return 0
	poa_args.ref = args_table['fasta']
	poa_args.out = args_table[(region,'poa')]
	poa_args.vcf = args_table[(region,'outvcf1')] + '.gz' ## longshot snvs vcf file
	poa_args.region = region
	#print(poa_args)

	print('\n\nrunning POA variant detection on haplotagged bam file',poa_args.bam,'in region',region,file=sys.stderr)
	varcaller = poacaller.VarCall(poa_args)
	varcaller.parse_region(region)
	varcaller.process_bam()
	varcaller.close_filehandles()
	poacaller.combinevcfs(poa_args.out,poa_args.ref,args_table['BCFTOOLS'])
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


def run_longshot(args_table,regions,filtered_args,PASS='first'):
	task_list = []
	logs_list = []
	for region in regions:
		print('\nrunning longshot',PASS,'pass for SNV calling',region,file=sys.stderr)

		if PASS=='second' and not os.path.isfile(args_table[(region,'haplobam')]): continue

		if PASS == 'first': longshot_cmd =args_table['LONGSHOT'] + ' ' + ' '.join(filtered_args) + ' --out_bam ' + args_table[(region,'haplobam')] + ' --out ' + args_table[(region,'outvcf1')] + ' --region ' + region
		elif PASS == 'second': longshot_cmd =args_table['LONGSHOT'] + ' ' + ' '.join(filtered_args) + ' -v ' + args_table[(region,'poavcf')]  + '.gz --out ' + args_table[(region,'outvcf2')] + ' --region ' + region + ' -D 100:500:50'
		#print(longshot_cmd)
		task_list.append(longshot_cmd)
		logs_list.append(args_table[(region,'log')])	

	process_tasks(task_list,logs_list,ncores=NCORES)
	task_list.clear()

def run_POA_python(args_table,regions):
	task_list = []
	logs_list = []
	for region in regions:
		if not os.path.isfile(args_table[(region,'haplobam')]): continue
		print('\nrunning POA variant calling on hap-reads',region,file=sys.stderr)
		poa_cmd ='python3.6 poacaller.py' + ' --bam ' + args_table[(region,'haplobam')] + ' --out ' + args_table[(region,'poa')] + ' --region ' + region + ' --vcf ' + args_table[(region,'outvcf1')] + '.gz' + ' --ref ' + args_table['fasta']
		task_list.append(poa_cmd)
		print(poa_cmd)
		logs_list.append(args_table[(region,'log')])	
	process_tasks(task_list,logs_list,ncores=NCORES)
	task_list.clear()



#################################################################################################################################	

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
		for contig in pyfasta.references: 
			rlength = pyfasta.get_reference_length(contig)
			regions.append(contig + ':1-' + str(min(rlength,250000)))
			#if len(regions) > 5: break
	else: regions.append(args_table['region'])

	counter=0
	for region in regions:
		args_table[(region,'haplobam')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.haplotag.bam'
		args_table[(region,'outvcf1')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.snvs.vcf'
		args_table[(region,'outvcf2')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.snvs_indels.vcf'
		args_table[(region,'poavcf')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.poa_filtered.vcf'
		args_table[(region,'poa')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.poa'
		args_table[(region,'log')] = args_table['OUTDIR'] + '/contig' + str(counter) + '.log'
		counter+=1

	##################################################################################################

	run_longshot(args_table,regions,filtered_args,PASS='first')
	print('gzipping and indexing vcf files',file=sys.stderr)
	for region in regions:
		subprocess.call('bgzip -f ' + args_table[(region,'outvcf1')],shell=True)
		subprocess.call('tabix -f ' + args_table[(region,'outvcf1')] + '.gz',shell=True)


	print('running POA consensus based variant detection',file=sys.stderr)
	run_POA_python(args_table,regions)

	print('filtering dense clusters of variants',file=sys.stderr)
	for region in regions:
		if not os.path.isfile(args_table[(region,'haplobam')]): continue
		filter_dense_variants(args_table[(region,'poa')]+'.final.diploid.vcf.gz',args_table[(region,'poavcf')])
		subprocess.call('bgzip -f ' + args_table[(region,'poavcf')],shell=True)
		subprocess.call('tabix -f ' + args_table[(region,'poavcf')] + '.gz',shell=True)
	
	run_longshot(args_table,regions,filtered_args,PASS='second')
	print('gzipping and indexing vcf files',file=sys.stderr)
	for region in regions:
		if not os.path.isfile(args_table[(region,'haplobam')]): 
			## if vcf file for a region is empty and no haplotagged bam file...
			subprocess.call('cp ' + args_table[(region,'outvcf1')] + '.gz'  + ' ' + args_table[(region,'outvcf2')] + '.gz',shell=True)
		else: 
			subprocess.call('bgzip -f ' + args_table[(region,'outvcf2')],shell=True)

		subprocess.call('tabix -f ' + args_table[(region,'outvcf2')] + '.gz',shell=True)
			
	print('combining individual vcf files into single output'); 
	concat_cmd = args_table['BCFTOOLS'] + ' concat -o ' + args_table['out'] + ' ' + ' '.join([args_table[(region,'outvcf2')] + '.gz' for region in regions])
	subprocess.call(concat_cmd,shell=True)
	
	##################################################################################################
	#futures=[]
	#with ThreadPoolExecutor(max_workers=NCORES) as executor:
	#	for region in regions: futures.append(executor.submit(call_POA,region,args_table))


if __name__ == "__main__":
	main()
