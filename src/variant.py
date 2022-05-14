import os,sys 
import copy
from .variant_functions import *
import statistics


class VCFrecord: ## standard VCF variant
	def __init__(self,chrom=0,pos=0,ref='.',alt='.',identifier='.', qual='.', filt='.', info='.',formatfield='GT',genotype=None):
		self.chrom=chrom
		self.pos = pos
		self.ref = ref
		self.alt=alt
		self.alts = [] ## list of alternate alleles
		self.identifier = identifier
		self.qual = qual
		if filt !='.': self.filt=[filt]
		else: self.filt= []
		self.info = info
		self.formatfield=formatfield
		self.genotype = genotype
		self.leftnormal = None ## left normalized representation of variant, namedtype Vartuple
		self.rightnormal = None ## right normalized representation of variant
		## phase information

	def normalize(self,refseq,offset):
		vt = Vartuple(pos=self.pos,ref=self.ref,alt=self.alt,alleles=[])
		self.leftnormal = normalize_left(vt,refseq,offset)
		#self.rightnormal = normalize_right(vt,refseq)
	
	def compare_to_Variant(self,other,refseq,offset,DEBUG=False,normalized=True): ## compare to 'Variant' record
		return self.compare(other.vcf,refseq,offset,normalized=False)

	def compare(self,other,refseq,offset,DEBUG=False):
		if self.leftnormal == None: self.normalize(refseq,offset)
		if other.leftnormal == None: other.normalize(refseq,offset)

		if self.leftnormal == other.leftnormal: return 1
		return compare_vartuples(self.leftnormal,other.leftnormal,refseq,offset=offset,normalized=True)
        	## determine rightmost possible position for each variant and check if there is overlap...

	def print(self,outfile=sys.stdout):
		self.genotype = '0/1'
		if len(self.filt) ==0: filt_string = '.'
		else: filt_string = ';'.join(self.filt)
		print(self.chrom,self.pos,self.identifier,self.ref,self.alt,self.qual,filt_string,self.info,'GT',self.genotype,sep='\t',file=outfile)


class Variant: ##variant extracted from POA consensus (aligned to reference using parasail or graph)
	def __init__(self,vtype=None,chrom=0,pos=0,ref='.',alt='.',pairwise=None,sindex=0,eindex=0,source='POA',haptag=None,maxcov=0):
		self.vtype = vtype
		self.chrom=chrom
		self.pos = pos; 

		## only for variant from parasail alignment
		self.pairwise = pairwise;
		self.start = sindex;  ## indices to pairwise vector
		self.end = eindex; 
		self.end_last = eindex
		self.ambig=0;

		self.mincov =1000; 
		self.mincov_ref = 1000;
		#self.mincov_other=1000
		self.maxcov = maxcov ## number of reads for haplotype
		self.ratio_alleles=0.0
		self.fraction=0.0
		self.product = 0.0
		self.flag = 'PASS'

		self.Wstart= 0; self.Wend=0; self.readcounts =None
		self.divergence = 0 ## divergence of consensus to reference
		self.variants_window =-1

		## not required 
		if self.pairwise != None:
			self.ref = pairwise[0][self.start:self.end];  ## same length as alt (includes '-' characters)
			self.alt = pairwise[2][self.start:self.end]; 
		else: 
			self.ref = ref
			self.alt = alt

		self.source = source
		self.validated = False
		self.genotype = None
		self.haptag = []
		if haptag != None: self.haptag.append(haptag) ## 0/1/2 or array
		self.matching_pair = None ## pointer to other variant that matches
		self.is_duplicate = False

		self.overlaps = [] ## list of other variants that overlap this variant 

		self.vcf = None

	def __lt__(self,other):
		if self.vcf != None and other.vcf != None:
			if self.vcf.pos == other.vcf.pos: return self.source < other.source
			return self.vcf.pos < other.vcf.pos
		else: 
			if self.pos == other.pos: return selfsource < other.source
			return self.pos < other.pos

	## also allows for identical bases in-between, ref and alt have to be equal length
	def splitMNP(self,DEBUG=True): ## return list of new variants or None | use deepcopy | ACG->TAA to A->T,C->A,G->A
		rl = len(self.ref)
		ra = len(self.alt)
		if rl ==1 or rl != ra or '-' in ''.join(self.ref) or '-' in ''.join(self.alt): return None
		nvars =0
		for i in range(rl):
			if self.ref[i] != self.alt[i]: nvars +=1
		if nvars== 1: return None

		if DEBUG: print('splitting MNP',self.ref,self.alt,nvars)
		varlist = []
		for i in range(rl):
			if self.ref[i] == self.alt[i]: continue
			newvar = copy.deepcopy(self) ## copy all attributes and change only three attributes
			newvar.pos = self.pos+i
			newvar.ref = self.ref[i:i+1]
			newvar.alt = self.alt[i:i+1]
			varlist.append(newvar)
		return varlist

	def compare(self,other,refseq,offset): ## both are Variant
		## determine if two variants are identical, useful for longshot comparison and merging haploid calls to diploid
		return self.vcf.compare(other.vcf,refseq,offset) ## use vcf class function

	## only function that depends on pairwise tuple
	def init_from_pairwise(self): ## pairwise is from Consensus (5-tuple)
		gaps = [0,0]
		s = -1; e=-1
		for i in range(self.start,self.end): 
			if self.pairwise[0][i] == '-': 
				gaps[0] +=1
				if s < 0: s = i
				e = i
			elif self.pairwise[2][i] == '-': 
				gaps[1] +=1
				if s < 0: s = i
				e = i
		
		if gaps[0] + gaps[1] == 0: 
			self.vtype = 'SNV'
		elif gaps[0] > 0 and gaps[1] ==0: ## insertion
			altseq = self.pairwise[2][s:e+1]
			self.vtype = 'INDEL+'
		elif gaps[0] == 0 and gaps[1] > 0: ## deletion
			altseq = self.pairwise[0][s:e+1]
			self.vtype = 'INDEL-'
		else: self.vtype ='COMPLEX'

		if 'INDEL' in self.vtype: 
			indel_length = gaps[0] + gaps[1]
			l=0; m = e+1
			right_limit = len(self.pairwise[0])
			while m < right_limit and altseq[l] == self.pairwise[0][m]: 
				l +=1 
				m +=1
				if l >= indel_length: l -= indel_length
			#print(altseq,s,gaps,'ambig',m-e-1)
			self.ambig = m-e-1
			self.end_last = m

		l = len(self.ref)
		if gaps[0] + gaps[1] ==0: 
			s = 0
			while s <l and  self.ref[s] == self.alt[s]: s +=1
			e = l-1
			while e >= s and  self.ref[e] == self.alt[e]: e -=1
		else:
			s=0 
			e = l-1
			while e >= s and  self.ref[e] == self.alt[e]: e -=1
		self.vcf = VCFrecord(chrom=self.chrom,pos=self.pos+s,ref = ''.join(self.ref[s:e+1]).replace('-',''),alt = ''.join(self.alt[s:e+1]).replace('-',''))

		for i in range(self.start,self.end_last):
			if i==self.end_last-1 and self.pairwise[3][i] == self.pairwise[4][i]: continue ## filter identical position 
			if self.pairwise[3][i] >= 0 and self.pairwise[3][i] < self.mincov_ref: self.mincov_ref = self.pairwise[3][i]
			if self.pairwise[4][i] >= 0 and self.pairwise[4][i] < self.mincov: self.mincov = self.pairwise[4][i]
			#if self.pairwise[5][i] >= 0 and self.pairwise[5][i] < self.mincov_other: self.mincov_other = self.pairwise[5][i]

	def calculate_statistics(self,nreads,minratio=4,mincount=3,outfile=sys.stdout):
		if self.mincov + self.mincov_ref > self.maxcov+1: print('coverage-BUG',self.mincov,self.mincov_ref,self.maxcov,file=outfile)
		elif self.mincov + self.mincov_ref < 1: 
			print('coverage-ZERO-BUG',self.mincov,self.mincov_ref,file=outfile)
			return 1
		self.ratio_alleles = float(self.mincov)/(self.mincov+self.mincov_ref)
		self.fraction = float(self.mincov)/nreads;
		self.product = self.fraction*self.ratio_alleles
		#self.ratio_alleles = float(self.mincov)/(max(0,self.mincov_ref)+0.1)

		## need to incorporate the fact that some variants are haptag '0', called using all reads from two haplotypes
		if self.haptag[0] == '0': ## diploid
			if self.fraction >= 0.15 and self.ratio_alleles >=0.3 and self.mincov >=mincount: self.flag='PASS'
			else: self.flag = 'filt'
		else: ## haploid
			if self.fraction >= 0.3 and self.ratio_alleles >=0.6 and self.mincov >=mincount: self.flag='PASS'
			else: self.flag = 'filt'

	def print_simpleformat(self,VCFformat=True,outfile=sys.stdout):
		print('var',self.pos,'amb',self.ambig,''.join(self.ref),''.join(self.alt),end='\t',file=outfile)
		if VCFformat: print('vcf:' + str(self.vcf.pos) + ':' + self.vcf.ref +':' + self.vcf.alt,end='\t',file=outfile)
		print('min',self.mincov_ref,self.mincov,self.flag,self.vtype,round(self.ratio_alleles,3),round(self.fraction,2),round(self.product,2),end=' ',file=outfile)
		if self.pairwise != None: 
			print('ref',self.pairwise[3][self.start:self.end_last],'alt',self.pairwise[4][self.start:self.end_last],end=' ',file=outfile)
			#print('alt1',self.pairwise[5][self.start:self.end_last],self.mincov_other,end=' ')
		print('hap',self.haptag,file=outfile) 

	## output cigar of full window, local cigar
	def print_VCFformat(self,REPLACE_AMBIG=True,outfile=sys.stdout,logfile=sys.stdout):
		if self.source == 'POA' and (self.mincov <3 or self.ratio_alleles < 0.5): return -1
		if self.vcf.ref == self.vcf.alt: 
			print('skipping vairant, identical ref and alt',self.vcf.pos,self.vcf.ref,file=logfile)
			return -1

		if len(self.vcf.filt) ==0: filt_string = '.'
		else: filt_string = ';'.join(self.vcf.filt)
		
		ref_list = list(self.vcf.ref)
		for i in range(len(ref_list)):
			if ref_list[i] != 'A' and ref_list[i] != 'C' and ref_list[i] != 'G' and ref_list[i] != 'T': 
				ref_list[i] = 'N'
				print('non-ACTG base in reference',self.vcf.ref,self.vcf.pos,file=logfile)
			else: pass
		refallele = ''.join(ref_list)	

		print(self.vcf.chrom,self.vcf.pos,self.vcf.identifier,refallele,self.vcf.alt,self.vcf.qual,filt_string,end='\t',sep='\t',file=outfile)

		if self.source == 'POA': 
			info = []
			if len(self.haptag) == 1 and self.haptag[0] == '1': self.vcf.genotype = '1|0'
			elif len(self.haptag) == 1 and self.haptag[0] == '2': self.vcf.genotype = '0|1'
			elif len(self.haptag) == 1 and self.haptag[0] == '0': 
				if self.ratio_alleles > 0.75: self.vcf.genotype= '1/1'
				else: self.vcf.genotype = '0/1'
			elif len(self.haptag) == 2: self.vcf.genotype = '1|1'
			else: self.vcf.genotype = '0/1'
			
			gtfields = ['GT']
			gtvalues = [self.vcf.genotype]
			gtfields.append('HP')
			gtvalues.append(self.haptag[0])
			if self.Wend > 0: 
				info.append(('W',str(self.Wstart)+'-'+str(self.Wend)))
				gtfields.append('WI')
				gtvalues.append(str(self.Wstart)+'-'+str(self.Wend))
			if self.readcounts != None: 
				#info.append(('RCOUNTS',str(self.readcounts[0]) + ','+str(self.readcounts[1]) + ','+str(self.readcounts[2])))
				gtfields.append('RC')
				if self.haptag[0] == '0': gtvalues.append(str(self.readcounts[0]+self.readcounts[1]+self.readcounts[2]))
				elif self.haptag[0] == '1': gtvalues.append(str(self.readcounts[1]))
				elif self.haptag[0] == '2': gtvalues.append(str(self.readcounts[2]))
				gtfields.append('RD')
				gtvalues.append(str(self.readcounts[0]+self.readcounts[1]+self.readcounts[2]))

			cigar='.'
			if self.vtype.split(':')[0] == 'SNV': 
				cigar = str(len(self.vcf.ref))+'M'
			elif self.vtype.split(':')[0] == 'INDEL':
				l = int(self.vtype.split(':')[1])
				if l < 0: cigar = str(-1*l) + 'D'
				else: cigar = str(l) + 'I'
			gtfields.append('RT') ## stats
			gtvalues.append(str(round(self.ratio_alleles,2)))
			gtfields.append('ST') ## stats
			gtvalues.append(str(round(self.fraction,2)) + ',' +str(round(self.divergence,3)) + ',' + str(self.variants_window)+','+cigar)

			print('SOURCE=POA_'+str(self.haptag[0]),':'.join(gtfields),':'.join(gtvalues),sep='\t',file=outfile)

			info += [('HP',self.ambig),('FL',self.flag),('VT',self.vtype),('RATIO',round(self.ratio_alleles,2))]
			info.append(('FRAC',round(self.fraction,2)))
			info.append(('PR',round(self.product,2)))
			info.append(('AC',str(self.mincov_ref)+',' +str(self.mincov)))
			info.append(('HAP',','.join(self.haptag)))
			info.append(('SOURCE',self.source))
			#if not self.validated: info.append(('INPUT','no'))
			#else: info.append(('INPUT','yes'))

			#print(gtfields,gtvalues)

			#print(';'.join([pair[0] + '=' + str(pair[1]) for pair in info]),'GT',self.vcf.genotype,sep='\t',file=outfile)

		else: ## from longshot or another caller

			gtfields = ['GT']
			gtvalues = [self.vcf.genotype.split(':')[0]]
			gtfields.append('HP')
			gtvalues.append('?')
			#gtfields.append('SRC')
			#gtvalues.append('INP')	
			print('SOURCE=INP',':'.join(gtfields),':'.join(gtvalues),sep='\t',file=outfile)

			info = ['VT=' + self.vtype,'W=' + str(self.Wstart)+'-'+str(self.Wend),'SOURCE=' + self.source]
			if len(self.vcf.info) > 1: info += self.vcf.info
			#print(';'.join(info),self.vcf.formatfield,self.vcf.genotype,sep='\t',file=outfile)
