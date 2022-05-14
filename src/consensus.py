import os,sys 
CIGAR_OP = ['M','I','D','N','S','H','P','=','X'];
import parasail
import statistics
from .variant import Variant,VCFrecord
from .variant_functions import *
from .spoawrapper import spoa_python,query_sequence_graph
from .graphbubble import CPath,Subgraph
from collections import Counter
from .alignment import *

def kmer_counts(sequence,n,klen=15,outfile=sys.stdout):
	if 'N' in sequence: return 100,[0,'NNN']
	kdict = {}
	count =0
	for i in range(n-klen):
		try: kdict[sequence[i:i+klen]] +=1
		except KeyError: kdict[sequence[i:i+klen]] =1
	for i in range(n-klen):
		print(kdict[sequence[i:i+klen]],end=',',file=outfile)
	print('',file=outfile)

## represent consensus information obtained from POA along with alignment
class Consensus:
	def __init__(self,ref_start=0,params=None,refseq=None):
		self.sequences=[]
		self.uniques = [] ## unique sequences for CCS reads
		self.refseq = refseq
		self.ref_offset =ref_start
		self.nreads =0
		self.haptag = -1 ## 0/1/2 
		self.otherCS = None
		self.varlist = [] ## list of variants
		self.varlist2= [] ## from bubbles
		self.subgraph = None ## Subgraph object that stores graph-based alignment/bubbles between POA consensus and reference
		self.params=params

		## objects returned from POA graph analysis
		self.consensus = None ## consensus and cvec have sam length
		self.cvec = None ## coverage of path corresponding to consensus, same length
		self.cvec_ref = None ## coverage of path corresponding to refseq, same length

		## alignment of consensus to information information
		self.perfectmatch = False
		self.alignment= None; ## returned object from parasail alignment
		self.cigartuples = None
		self.divergence =0.0;
		self.filter = False;
		self.pairwise = None ## three strings with full alignment
		self.ref_aln_start = ref_start ## position on reference genome where the consensus aligns to (start position)
		self.ref_aln_end = ref_start  ## position where alignment ends


	def _print_(self,pflag='',outfile=sys.stdout):
		if self.alignment == None: 
			if pflag != 'all0': print(pflag,self.nreads,'.',file=outfile)
		elif self.perfectmatch: print(pflag,self.nreads,'perfect',file=outfile)
		else: 
			largecigars=0
			for seq in self.sequences:
				for op in seq[-1]: 
					l = len(op)
					if (op [l-1] == 'D' or op[l-1] == 'I') and int(op[0:l-1]) >= 50: largecigars +=1

			print('\n',pflag,self.nreads,self.ref_aln_start,self.ref_aln_end+1,self.ref_aln_end+1-self.ref_aln_start,end=' ',file=outfile)
			print(':'.join(self.cigartuples),'divergence',round(self.divergence,3),largecigars,file=outfile)

			if self.divergence >=0.03 or largecigars >=3:
				for seq in self.sequences: print ('high',seq[-1],len(seq[1]),file=outfile)
				#kmer_counts(self.refseq,len(self.refseq))
				#print(self.consensus)
			
			if self.pairwise == None: print('pairwise vector missing',file=outfile)
			else:	print(''.join(self.pairwise[0]) + '\n' + ''.join(self.pairwise[1]),file=outfile) ## ref vs alt

			#for v in self.varlist: 
			#	v.print_simpleformat()  ## print variants
			#print('')

			for b in self.subgraph.bubble_list: 
				#print('leftnormal',b.var.vcf.leftnormal)
				b.print_bubble(outfile=outfile)

	## parse cigar -> pairwise -> get variants
	def extract_variants(self,minsep=1,SPLIT_MNP=False,chrom=0,min_padding=10):  
	## pairwise[0] is reference, ## pairwise[1] is edit, pairwise[3] is ref_coverage, [4] -> alt_coverage
		if len(self.cigartuples) ==1 or self.alignment==None: return []
		self.pairwise = self.cigar2pairwise()
		flag=0; 
		start = 0; end=0
		shift_left_ref=0; ## keep track of how many insertions in reference
		previous_boundary = 0
		variants =0
		varlist = []
		l =len(self.pairwise[0])
		for i in range(l):
			if flag ==0 and self.pairwise[1][i] != '.': 
				start =i-1; flag = 1
				#print('varpos',start,self.ref_start,shift_left_ref)
				varpos = start+self.ref_aln_start+1-shift_left_ref
			elif flag > 0 and self.pairwise[1][i] == '.':  ## if two SNVs next to each other they will be merged
				if flag < minsep: flag +=1
				else: 
					end=i+1
					var =Variant(chrom=chrom,pos=varpos,pairwise=self.pairwise,sindex=start,eindex=end,source='POA',haptag=self.haptag,maxcov=self.nreads); 
					var.prevb = previous_boundary
					if variants > 0: varlist[variants-1].nextb = start-1 

					if SPLIT_MNP: tempvarlist = var.splitMNP()
					else: tempvarlist = None
					if tempvarlist != None and i >= min_padding and i < l-min_padding:
						varlist += tempvarlist ## split variant and add multiple variants to list
						variants += len(tempvarlist)
					elif i >= min_padding and i < l-min_padding:	
						varlist.append(var);
						variants +=1
					flag = 0
					previous_boundary = end
			if self.pairwise[0][i] == '-': shift_left_ref +=1
		if variants > 0: varlist[variants-1].nextb = i
		return varlist
			

	def cigar2pairwise(self):
		if self.alignment == None: return [-1,[]];
		ref = []; comparison = []; query = []; 
		coverage = []; 
		coverage_other = []; 
		coverage_ref = []
		#print('debug',self.cvec,self.cvec_ref,len(self.refseq),len(self.cvec_ref))
		length_p=0; ## pairwise alignment vectors
		cigar = (self.alignment).cigar;
		lr=cigar.beg_ref; lq=cigar.beg_query;
		self.ref_aln_start += lr ## alignment starts at position 'lr' relative to 'refseq' instead of 0
		n= len(cigar.seq)
		for i in range(n):
			op = CIGAR_OP[cigar.seq[i] &15]; 
			ol = cigar.seq[i]>>4; 
			if op =='D' and (i ==0 or i == n-1): op = 'S'; ## treat it like soft-clip

			if op == 'S': 
				if i == n-1: break; # dont increment lr 
				lr += ol;
				self.ref_aln_start += ol # since we don't add columns to pairwise
			elif op == 'D': 
				for j in range(ol): 
					ref.append(self.refseq[lr+j]); comparison.append('D');  query.append('-')
					coverage.append(-1); coverage_ref.append(self.cvec_ref[lr+j])
					coverage_other.append(-1)
					length_p +=1;
				lr += ol; 
			elif op == 'I':
				## skip insertion at first or last position in alignment 
				if i > 0 and i < n-1: 
					for j in range(ol): 
						ref.append('-'); query.append(self.consensus[lq+j]); comparison.append('I')
						coverage.append(self.cvec[lq+j]); coverage_ref.append(-1)
						length_p +=1;
				lq += ol; 
			elif op == '=': 
				for j in range(ol): 
					ref.append(self.refseq[lr+j]); comparison.append('.'); query.append(self.refseq[lr+j])
					coverage.append(self.cvec[lq+j]); coverage_ref.append(self.cvec_ref[lr+j]);
					length_p +=1;
				lr += ol; lq += ol;
			elif op == 'X': 
				for j in range(ol):
					ref.append(self.refseq[lr+j]); query.append(self.consensus[lq+j]); comparison.append('X')
					coverage.append(self.cvec[lq+j]); coverage_ref.append(self.cvec_ref[lr+j]);
					length_p +=1;
				lr +=ol; lq +=ol;	
		self.ref_aln_end += lr
		pairwise = [ref,comparison,query,coverage_ref,coverage,coverage_other]
		#for i in range(len(pairwise[0])): print(pairwise[0][i],pairwise[3][i],pairwise[2][i],pairwise[4][i],pairwise[1][i],i)
		return pairwise
		#print('parasail',','.join(self.cigartuples),'BR',cigar.beg_ref,'BQ',cigar.beg_query,self.alignment.score,self.ref_aln_start,self.ref_aln_end);

	def callPOA(self,haptag='0',chrom=0,mincov=3,compareOther=False,use_alignment=False,SPLIT_MNP=False,normalize=False,outfile=sys.stdout): 
		## if normalize, use own left normalization
		self.haptag = haptag
		#for seq in self.sequences: 
		#	if seq[1] == '': print('BUG...empty sequence',seq)

		print ("SPOA on reads for haptag",haptag,self.nreads,file=outfile)
		consensus,cvec,self.cvec_ref,paths = self.spoa_prepare(haptag=haptag,outfile=outfile);
		if cvec == None: return 0;
		#print('cons',paths[0],'\n',paths[1],'\n',consensus,'\n',cvec) ## DEBUGverbo

		## trim low-coverage ends of consensus before calling alignment-based variant detection
		nc = len(consensus);
		thresh = max(mincov,int(0.4*self.nreads));
		s = 0; e = nc-2;
		while s < e and cvec[s] < thresh: s +=1;
		while e > s and cvec[e] < thresh: e -=1;
		self.consensus=consensus[s:e+1]
		self.cvec = cvec[s:e+1]

		nvars = -1
		if self.refseq.find(self.consensus) >=0: 
			self.perfectmatch =True
			return 0
		else: 
			self.alignment,self.cigartuples,self.divergence = align_parasail(self.consensus,self.refseq)
			if self.divergence >= 0.25 or self.alignment.score <=0: self.filter=True
			varlist1 = self.extract_variants(chrom=chrom,SPLIT_MNP=SPLIT_MNP,min_padding=10) ## option 1
			for var in varlist1: 
				var.init_from_pairwise() ## initialize variant values using pairwise tuple
				#var.calculate_statistics(self.nreads)
			nvars = len(varlist1)

		## paths = (cons_path,ref_path), graph-based variant detection
		self.subgraph=Subgraph(paths[0],paths[1],outfile=outfile)
		#if self.subgraph.trimmed ==1:
		#for seq in self.sequences: print('debug',seq[0],seq[-1])
		self.subgraph.extract_bubbles_graph(self.ref_offset,self.nreads) 
		varlist2 = []
		for bubble in self.subgraph.bubble_list:
			var = bubble.bubble_to_variant(chrom=chrom,haptag=haptag,maxcov=self.nreads)
			if SPLIT_MNP: multivars = var.splitMNP()
			else: multivars = None
			if multivars == None: varlist2.append(var)
			else: varlist2 += multivars

		lvars2 = len(varlist2)
		for var in varlist2:
			var.calculate_statistics(self.nreads,outfile=outfile)
			var.divergence=self.divergence
			var.variants_window = lvars2
			vt = Vartuple(pos=var.pos,ref=var.ref,alt=var.alt,alleles=[])
			if normalize: 
				leftnormal = normalize_left(vt,self.refseq,self.ref_offset)
				var.vcf = VCFrecord(chrom=chrom,pos=leftnormal.pos,ref =leftnormal.ref,alt =leftnormal.alt)
			else:
				var.vcf=VCFrecord(chrom=chrom,pos=var.pos,ref=var.ref,alt=var.alt)
			print ("bubble-derived",end=' ',file=outfile); var.print_simpleformat(outfile=outfile)

		#for var in varlist1: print ("alignment-derived",end=' '); var.print_simpleformat()
		if use_alignment: self.varlist =varlist1
		else: self.varlist =varlist2

		if nvars >= 0 and nvars != self.subgraph.nbubbles: print('difference in methods',nvars,self.subgraph.nbubbles,file=outfile)


	## ADD_FIRST, CIGAR_FIRST,parameters
	def spoa_prepare(self,haptag,DEBUG=True,outfile=sys.stdout):
		if haptag =='0': diploid = True ## joint POA
		else: diploid = False

		nu = len(self.uniques)
		if diploid == False and nu > 0 and self.uniques[0][0] >= 3 and self.uniques[0][0]*2 > self.nreads:
			if self.refseq.find(self.uniques[0][1]) >=0: 
				#print('using refseq for spoa perfect',nu,self.uniques[0])
				if DEBUG: print('perfect hit without POA',self.uniques[0][0],self.nreads,file=outfile)
				return self.uniques[0][1],None,None,None

		if nu <1 or self.nreads-nu <=2: ## only small difference between total_reads and unique_reads 
			return spoa_python(self.sequences,weighted=False,refseq=self.refseq,ADD_FIRST=True,DIPLOID=diploid,outfile=outfile);
			print('using refseq for spoa perfect',nu)
			## consensus,coverage_consensus,coverage_refseq
		else:
			spoa_items= spoa_python(self.uniques,weighted=True,refseq=self.refseq,ADD_FIRST=False,DIPLOID=diploid,outfile=outfile);
			return spoa_items
			print('spoa_consensus',len(spoa_items[0]),file=outfile)
			print('saved ops',self.nreads-nu,len(self.refseq),file=outfile);

	def get_coverage_otherhap(self,other,POA): ## add sequences from 2nd haplotype to graph for this haplotype and get coverage for consensus vector
		nu = len(other.uniques)
		if nu < 1 or other.nreads-nu <= 2: 
			for i in range(other.nreads): POA[0].add_sequence(POA[1],other.sequences[i][1],1);
		else: 
			for read in other.uniques:  POA[0].add_sequence(POA[1],read[0],read[1])
		cvec =  [int(edge.weight()/2) for edge in POA[2]]; cvec.append(0)
		return cvec

	def get_uniques(self):
		self.uniques = []
		table = {}
		for read in self.sequences:
			try: table[read[1]][0] +=1 
			except KeyError: table[read[1]] = [1,read[-1]] # count,cigar tuple
		for seq in table.keys(): self.uniques.append([table[seq][0],seq,table[seq][1]]) ## (count,sequence,cigar)
		self.uniques.sort(reverse=True)

		#self.uniques = list(Counter([read[1] for read in self.sequences]).items()) ## get unique sequence counts
