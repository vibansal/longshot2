import os,sys 
CIGAR_OP = ['M','I','D','N','S','H','P','=','X'];
import parasail

## parasail starts with gapO for a gap of length 1, instead of gapO+gapE, 
##so it is important to use higher penalty since other methods use a+ b.k where k is gap length
def align_parasail(consensus,refseq,match=2,mismatch=-4,gapO=6,gapE=2): 
	aln_matrix = parasail.matrix_create("ACGT", 2, -4);
	alignment=parasail.sw_trace(consensus,refseq,6,2,aln_matrix);
	#alignment=parasail.sg_dx_trace(consensus,refseq,6,3,aln_matrix); #do not penalize gaps at begin and end of s1/query
	cigar=alignment.cigar; Q = alignment.query;
	cigartuples = [];
	n=len(cigar.seq);
	edit=0; length=0;
	for i in range(n):
		op = CIGAR_OP[cigar.seq[i] &15];
		ol = cigar.seq[i]>>4; 
		#print('cigar',ol,op)
		if op == 'X' or op == 'D' or op == 'I': edit+= 1; ## single event
		if op =='D' and (i ==0 or i == n-1): edit -=1
		if op == '=': length += ol
		if op != 'S': cigartuples.append(str(ol)+op);
	divergence=float(edit)/(edit+length);
	#if divergence >= 0.25 or edit < 1 or alignment.score <=0: self.filter=True;
	return alignment,cigartuples,divergence
	
def cigar2pairwise_1(alignment,refseq,cvec_ref,cvec):
	if alignment == None: return [-1,[]];
	ref = []; comparison = []; query = []; 
	coverage = []; 
	coverage_other = []; 
	coverage_ref = []
	length_p=0; ## pairwise alignment vectors
	cigar = (alignment).cigar;
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
			self.ref_aln_start += lr # since we don't add columns to pairwise
		elif op == 'D': 
			for j in range(ol): 
				ref.append(refseq[lr+j]); comparison.append('D');  query.append('-')
				coverage.append(-1); coverage_ref.append(cvec_ref[lr+j])
				coverage_other.append(-1)
				length_p +=1;
			lr += ol; 
		elif op == 'I':
			## skip insertion at first or last position in alignment 
			if i > 0 and i < n-1: 
				for j in range(ol): 
					ref.append('-'); query.append(consensus[lq+j]); comparison.append('I')
					coverage.append(cvec[lq+j]); coverage_ref.append(-1)
					length_p +=1;
			lq += ol; 
		elif op == '=': 
			for j in range(ol): 
				ref.append(refseq[lr+j]); comparison.append('.'); query.append(refseq[lr+j])
				coverage.append(cvec[lq+j]); coverage_ref.append(cvec_ref[lr+j]);
				length_p +=1;
			lr += ol; lq += ol;
		elif op == 'X': 
			for j in range(ol):
				ref.append(refseq[lr+j]); query.append(consensus[lq+j]); comparison.append('X')
				coverage.append(cvec[lq+j]); coverage_ref.append(cvec_ref[lr+j]);
				length_p +=1;
			lr +=ol; lq +=ol;	
	self.ref_aln_end += lr
	pairwise = [ref,comparison,query,coverage_ref,coverage,coverage_other]
	#for i in range(len(pairwise[0])): print(pairwise[0][i],pairwise[3][i],pairwise[2][i],pairwise[4][i],pairwise[1][i],i)
	return pairwise
	#print('parasail',','.join(self.cigartuples),'BR',cigar.beg_ref,'BQ',cigar.beg_query,self.alignment.score,self.ref_aln_start,self.ref_aln_end);
