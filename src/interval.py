from __future__ import print_function
import os,sys
from .consensus import Consensus
from .variant import Variant,VCFrecord
import math

def kmer_complexity(sequence,n,klen=3,pflag=False):
	if 'N' in sequence: return 100,[0,'NNN']
	kdict = {}
	count =0
	for i in range(n-klen):
		try: kdict[sequence[i:i+klen]] +=1
		except KeyError: kdict[sequence[i:i+klen]] =1
		count +=1
	entropy = 0.0
	maxcount = [0,None]
	for kmer in kdict.keys(): 
		c = kdict[kmer]
		#if c < 2: continue
		f = float(c)/count
		entropy += f*math.log(f)
		if c > maxcount[0]: maxcount = [c,kmer]
	if pflag: print('ent',round(-1*entropy,2),maxcount)
	return -1*entropy,maxcount

class Interval(object):
	def __init__(self,chrom,left,right,params=None,seq=None):
		self.chrom=chrom
		self.left = left
		self.right =right 
		self.seq = seq ## reference
		self.leftpad = ""
		self.rightpad = ""
		self.paddedseq = None # with padding for variant analysis
		self.paddedstart = 0 ## start of paddedseq in reference genome
		self.params=params
		
		self.hapcounts = [0,0,0] ## count of reads assigned to each haplotype (?/1/2)
		self.nreads=0;
		self.partial =0
		self.MQcounts = [0,0,0,0] ## <10,10-19,20-29,30+

		self.CS = None; self.CS1 = None; self.CS2 = None

		self.POAoutcomes = [None,None,None] ## whether POA was called and if it yielded any variant(s) 
		self.varlist = [] ## merged list of variants from hap1,hap2,all

		self.input_variants = []

	def filter_outlier_reads(self):
		prev=0; 
		mean=0.0
		for read in self.CS.sequences: 
			c = read[5][0]+read[5][1];
			if prev > 0: mean += c-prev
			prev=c
		mean /= self.nreads-1
		print('mean-sep',mean)	
		prev=0; 
		for i in range(1,self.nreads): 
			read = self.CS.sequences[i]
			c = read[5][0]+read[5][1];
			print('read',round(mean,2),read[2],round(read[3][0],3),round(read[3][1],3),round(read[3][2],3),c-prev,read[5]) #read.subcigar);
			if c-prev > 6*mean and prev > 0: print('flag',c-prev,mean,end='\n') 
			prev=c
	
	def combine_variants(self,outfile=sys.stdout,normalize=False): ## don't ignore combined set 
		#print('merging variants from two haplotypes',v1,v2,v0)
		self.varlist = self.CS1.varlist + self.CS2.varlist + self.CS.varlist
		self.varlist.sort()

		for var in self.varlist: 
			if normalize: var.vcf.normalize(self.paddedseq,offset=self.paddedstart) 
			var.Wstart = self.left
			var.Wend = self.right
			var.readcounts = self.hapcounts

	def populate_slices(self,readlist,refseq,outfile=sys.stdout):

		self.CS = Consensus(self.left,params=self.params,refseq=self.seq)
		self.CS1 = Consensus(self.left,params=self.params,refseq=self.seq)
		self.CS2 = Consensus(self.left,params=self.params,refseq=self.seq)
		for read in readlist:
			read.getslice(self.left,self.right,min_overlap=-1,outfile=outfile);
			if read.partial: self.partial +=1
			if read.subread == None: continue;
			mqbin = min(3,int(read.MQ/10)); self.MQcounts[mqbin] +=1
			if read.MQ < self.params.minMQ: continue
			if self.params.filterreads and refseq.find(read.subread) >= 0: continue; ## filter exact matches for Illumina 
			self.CS.sequences.append([read.readid,self.leftpad+read.subread+self.rightpad,read.haptag,read.estats,read.strand,read.indelcount,read.partial,read.subcigar]); 
			self.nreads +=1;

		self.CS.sequences.sort(key= lambda x:x[5][0]+x[5][0]); ## indel counts of read within window, CLR reads
		#self.filter_outlier_reads()
		## partial reads should be last in the sorted list, not implemented

		for i in range(self.nreads):
			self.CS.nreads +=1
			if self.CS.sequences[i][2]==0: self.hapcounts[0] +=1; 
			elif self.CS.sequences[i][2]==1: self.CS1.sequences.append(self.CS.sequences[i]); self.hapcounts[1] +=1;
			elif self.CS.sequences[i][2]==2: self.CS2.sequences.append(self.CS.sequences[i]); self.hapcounts[2] +=1;
		self.CS1.nreads = self.hapcounts[1]
		self.CS2.nreads = self.hapcounts[2]

	def find_variants_POA(self,DEBUG=False,outfile=sys.stdout):
		minreads=self.params.minreads
		nc = self.seq.count('N'); 
		if nc > 0 or self.nreads > self.params.maxdepth or self.nreads < minreads: 
			if DEBUG: print('did not call POA.... reads',self.nreads,'#N',nc,self.nreads,self.left,self.right,file=outfile);
			return None

	        #if DEBUG: print ("\n########  SPOA on reads for haplotype2",self.hapcounts[2],'unassigned',self.hapcounts[0],self.left,self.right)
		n = self.nreads
		if (float(self.hapcounts[0])/n > 0.75 or self.hapcounts[1] < minreads or self.hapcounts[2] < minreads) or self.params.nohaps:
			if n >= minreads: self.CS.callPOA(haptag='0',chrom=self.chrom,outfile=outfile)
			else: 
				if DEBUG: print('too few sequences for POA',file=outfile);
		else:
			if float(self.hapcounts[0])/n > 0.25 and self.hapcounts[0] >= minreads: self.CS.callPOA(haptag='0',chrom=self.chrom,outfile=outfile)

			if self.hapcounts[1] >= minreads: 
				#self.CS2.otherCS = self.CS1
				if self.params.cluster: self.CS1.get_uniques()
			if self.hapcounts[2] >= minreads: 
				#self.CS1.otherCS = self.CS2
				if self.params.cluster: self.CS2.get_uniques()

			if self.hapcounts[1] >= minreads: self.CS1.callPOA(haptag='1',chrom=self.chrom,compareOther=self.params.compare,outfile=outfile)

			if self.hapcounts[2] >= minreads: self.CS2.callPOA(haptag='2',chrom=self.chrom,compareOther=self.params.compare,outfile=outfile)

	def _print_(self,nvar,outfile=sys.stdout):
		flag=0
		if (self.CS.alignment != None and not self.CS.perfectmatch): flag +=1
		if (self.CS1.alignment != None and not self.CS1.perfectmatch): flag +=1
		if (self.CS2.alignment != None and not self.CS2.perfectmatch): flag +=1
		#	for row in pybed.fetch(self.chrom,self.left,self.right): print('giab',str(row))
		print(flag,'window',self.left,self.right,self.nreads,'hap-counts',self.hapcounts,end=' ',file=outfile)
		print(self.seq[0:20],'...',self.seq[len(self.seq)-20:],'pt',self.partial,','.join([str(c) for c in self.MQcounts]),file=outfile);
		if flag > 0 or nvar > 0:
			self.CS._print_(pflag='all0',outfile=outfile); 
			self.CS1._print_(pflag='hap1',outfile=outfile); 
			self.CS2._print_(pflag='hap2',outfile=outfile); 
		#else: print('novar')


