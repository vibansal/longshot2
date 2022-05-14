import pysam,os,sys,glob,subprocess
import math,argparse,time,random

def read_error_stats(read,pflag=True,minmq=0,min_length=4000):
	read_length = read.query_alignment_length;
	if read.flag > 1024 or read.mapping_quality < minmq or read_length < min_length: return [0.01,0.01,0.01,read_length];
	M =0.0; D=0.0; I=0.0; S =0; D2 =0.0; I2=0.0;
	Dlist = [0,0,0,0,0]; Ilist = [0,0,0,0,0]; 
	mismatches = 0.0; insertions =0.0; deletions=0.0;
	for cigar in read.cigartuples: 
		if cigar[0] == 8: M += cigar[1];
		if cigar[0] == 1: 
			if cigar[1] ==1: I +=1;  
			if cigar[1] < 5: Ilist[cigar[1]] +=1;
		if cigar[0] == 2: 
			if cigar[1] ==1: D +=1;  
			if cigar[1] < 5: Dlist[cigar[1]] +=1;
		if cigar[0] == 4: S += cigar[1];

	ref_length = read.reference_end-read.reference_start;
	excess_ins = float(read_length-ref_length)/read_length;
	try: 
		RG = read.get_tag('RG');	
	except KeyError: RG = None;
	if pflag:
		print (round(M/read_length,4),round(D/read_length,4),round(I/read_length,4),round(excess_ins,4),end=' ');
		#if RG != None: print RG,
		print (read_length,ref_length,'soft-clipped',S,end=' ')
		print ('subs',M,'del',D,'ins',I,Dlist[1:],Ilist[1:],read.query_name);
	return [M/read_length,D/read_length,I/read_length,read_length];	

