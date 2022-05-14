import os,sys
import collections
import copy

Vartuple = collections.namedtuple('Vartuple','pos ref alt alleles')   
## (43444,ACC,A,[]) 
## alleles is non-empty for multi-allelic variants and alt = '.', alleles[0] = first_alt_allele, alleles[1] = second_alt_allele...

"""
functions for normalizing vcf bi-allelic variants, input is Vartuple
functions don't check for errors in variant ref allele
variant position is using +1 offset rather than 0-offset
"""

## returns vartuple
def normalize_right(vcfvar,refseq,offset=0,DEBUG=False): ## reverse algorithm, delete same base from left if identical, pad one base to right 
	rlimit=len(refseq)
	lr = len(vcfvar.ref); la = len(vcfvar.alt)
	cref = list(vcfvar.ref); calt = list(vcfvar.alt)
	lr_fixed = lr
	changes =1
	shiftright = 0
	while changes > 0:
		changes=0
		if lr > 0 and la > 0 and cref[0] == calt[0]:  ## delete same base from left if identical
			cref = cref[1:]
			calt = calt[1:]
			changes +=1
			lr -=1; la -=1
		if lr ==0 or la == 0: ## empty allele, pad one base to right from reference  TCCC[pad] T
			if DEBUG:print('c1',cref,calt,shiftright,vcfvar.pos,'offset',offset,lr_fixed)
			if vcfvar.pos-offset+lr_fixed+shiftright >= rlimit: 
				print('ERROR normalize_right',vcfvar,refseq,offset,lr_fixed,shift_right)
				return vcfvar
			nextbase = refseq[vcfvar.pos-offset+lr_fixed+shiftright]
			cref.append(nextbase)
			calt.append(nextbase)
			lr +=1; la +=1
			shiftright +=1
	shiftleft =0
	while lr >= 2 and la >= 2 and cref[lr-1-shiftleft] == calt[la-1-shiftleft]:
		lr -=1; la -=1 
		shiftleft +=1

	if DEBUG: print('final',vcfvar.pos,vcfvar.pos+shiftright,cref[0:lr-shiftleft],calt[0:la-shiftleft],shiftright,shiftleft,refseq)
	normvar = Vartuple(pos=vcfvar.pos+shiftright,ref=''.join(cref[0:lr-shiftleft]),alt=''.join(calt[0:la-shiftleft]),alleles=[])
	return normvar#,shiftright,shiftleft


## left justify and minimal length, algorithm from https://academic.oup.com/bioinformatics/article/31/13/2202/196142
def normalize_left(vcfvar,refseq,offset=0,DEBUG=False): ## vcfvar.pos-offset is 1-indexed into refseq
	#print('test',vcfvar.pos,vcfvar.ref,vcfvar.alt,refseq[vcfvar.pos-offset-1:vcfvar.pos-offset+10])
	lr = len(vcfvar.ref); la = len(vcfvar.alt)
	#if lr == 1 and lr == la: return copy.deepcopy(vcfvar)
	cref = list(vcfvar.ref); calt = list(vcfvar.alt) 
	changes =1
	shiftleft = 0
	while changes > 0:
		changes=0
		if lr > 0 and la > 0 and cref[lr-1] == calt[la-1]:  ## delete same base from right if identical
			cref.pop()
			calt.pop()
			changes +=1
			lr -=1; la -=1
		if lr ==0 or la == 0: ## empty allele, pad one base to left
			prevbase = refseq[vcfvar.pos-offset-2-shiftleft]
			cref = [prevbase] + cref
			calt = [prevbase] + calt
			#print(vcfvar.pos-1-shiftleft,prevbase,'c1',cref,calt)
			shiftleft +=1
			changes +=1
			lr +=1; la +=1
		if DEBUG: print('strings',cref,calt,changes,lr,la)
			 
	shiftright =0
	while lr >= 2 and la >= 2 and cref[shiftright] == calt[shiftright]:
		lr -=1; la -=1 
		shiftright +=1

	normvar = Vartuple(pos=vcfvar.pos-shiftleft+shiftright,ref=''.join(cref[shiftright:]),alt=''.join(calt[shiftright:]),alleles=[])
	if DEBUG: print('final',normvar,shiftleft,shiftright)
	return normvar#,shiftleft,shiftright

def compare_vartuples(var1,var2,refseq,offset=0,normalized=False,DEBUG=False,merge=False,maxlength=60):
	if len(var1.ref) >= 60 or len(var2.ref) >= maxlength: 
		print('very long variants.. skipping')
		return 0
	if normalized:
		var1_ln = var1
		var2_ln = var2
	else:
		var1_ln = normalize_left(var1,refseq,offset=offset)
		var2_ln = normalize_left(var2,refseq,offset=offset)
		#print(var1_ln,var2_ln)
	var1_endpos = var1_ln.pos + len(var1_ln.ref)-1
	var2_endpos = var2_ln.pos + len(var2_ln.ref)-1		
	if DEBUG: print(var1_ln,var2_ln,offset,'compare_var func')

	if var1_ln == var2_ln: 
		print('identical variants')
		return 1
	elif var2_ln.pos >= var1_ln.pos and var2_endpos <= var1_endpos:
		delta = var2_ln.pos-var1_ln.pos
		if var2_ln.ref == var1_ln.ref[delta:delta+len(var2_ln.ref)] and var2_ln.alt == var1_ln.alt[delta:delta+len(var2_ln.alt)]:
			print('var1 contains var2',var1_ln,var2_ln)
			return 4
	elif var1_ln.pos >= var2_ln.pos and var1_endpos <= var2_endpos:
		delta = var1_ln.pos-var2_ln.pos
		if var1_ln.ref == var2_ln.ref[delta:delta+len(var1_ln.ref)] and var1_ln.alt == var2_ln.alt[delta:delta+len(var1_ln.alt)]:
			print('var2 contains var1',var1_ln,var2_ln)
			return 4

	if (var2_ln.pos <= var1_ln.pos and var2_endpos >= var1_ln.pos) or (var1_ln.pos <= var2_ln.pos and var1_endpos >= var2_ln.pos):
		print('overlapping-simple',var1_ln,var2_ln)
		return 4
	else: ## last resort
		var1_rn = normalize_right(var1,refseq,offset=offset)
		var2_rn = normalize_right(var2,refseq,offset=offset)
		poslist = [(var1_ln.pos,'l',1),(var1_rn.pos+len(var1_rn.ref)-1,'r',1),(var2_ln.pos,'l',2),(var2_rn.pos+len(var2_rn.ref)-1,'r',2)]
		poslist.sort()
		if poslist[0][2] != poslist[1][2]: 
			print('overlapping-complex',var1_ln,var2_ln,poslist)
			return 4
		else: 
			#print('distinct',poslist)
			return 0

## combining variants on same haplotype to return a single bi-allelic variant
def combine_vartuples(varlist,refseq,offset=0,DEBUG=False): 
	pass;

## returns new vartuple with 'alleles' set to non-empty
## this function should be called after compare_vartuples, assume that there is overlap
## first align the position of the variants, then the length of the reference allele
## this is for merging variants from different callers
def merge_vartuples(varlist,refseq,offset=0,DEBUG=False): ## can merge an arbitrary number of bi-allelic variants
	if DEBUG: print(type(varlist[0]),file=sys.stderr)	
	reflist = []
	lengths = []
	nvars = len(varlist)
	minpos=varlist[0].pos
	for i in range(nvars):
		reflist.append(list(varlist[i].ref)); 
		lengths.append(len(varlist[i].ref))
		if varlist[i].pos < minpos: minpos = varlist[i].pos
	maxref = [0,0]

	## pad all variants to have same start position
	for i in range(nvars):
		l = lengths[i] + varlist[i].pos-minpos 
		if l > maxref[0]: maxref[0] =l; maxref[1] = i
		lengths[i] =l
		if varlist[i].pos == minpos: continue
		padding = refseq[minpos-offset-1:varlist[i].pos-offset-1]
		reflist[i] = list(padding) + reflist[i]

	## all variants have same length reference allele, alt-allele padded appropriately
	alleles = [] # ''.join(reflist[maxref[1]])] ## reference allele
	for i in range(nvars):
		alleles.append(refseq[minpos-offset-1:varlist[i].pos-offset-1] + varlist[i].alt + ''.join(reflist[maxref[1]][lengths[i]:maxref[0]] ))

	mergedvar = Vartuple(pos=minpos,ref=''.join(reflist[maxref[1]]),alt='.',alleles=alleles)
	if DEBUG: print('vars',varlist,'\n','final',mergedvar.pos,mergedvar.ref,mergedvar.alleles,file=sys.stderr)
	return mergedvar
	

if __name__ == "__main__":
	refseq = 'GCTTATCACACACACACACACACACCAAAAAGTG'
	## chr20   102609  .       CTGTGTGTGTGTGTGTGTGTGTG C, 102609:CTG:C
	## 157927  .       ATATATGTG       A 157929  .       ATATGTGTG       A
	## 374905:A:AC and 374906:A:C
	varlist= []
	varlist.append(Vartuple(pos=10,ref='CTGTGTGTGTGTGTGTGTGTGTG',alt='C',alleles=[])); varlist.append(Vartuple(pos=10,ref='CTG',alt='C',alleles=[]))
	#varlist.append(Vartuple(pos=30,ref='A',alt='AC',alleles=[])); varlist.append(Vartuple(pos=31,ref='A',alt='C',alleles=[]))
	#varlist.append(Vartuple(pos=5,ref='TCACACA',alt='T',alleles=[])); varlist.append(Vartuple(pos=11,ref='A',alt='C',alleles=[])); varlist.append(Vartuple(pos=8,ref='CA',alt='CACACA',alleles=[]))
	merge_vartuples(varlist,refseq,DEBUG=True)	

	sys.exit()

	# execute only if run as a script
	#var = Vartuple(pos=11,ref='ACA',alt='CGT'); var1 = Vartuple(pos=13,ref='A',alt='T')
	var = Vartuple(pos=11,ref='ACA',alt='A'); 
	var1 = Vartuple(pos=15,ref='AC',alt='A',alleles=[])
	compare_variants(var,var1,refseq)
	#print(refseq[10:])
	#normalize_left(var,refseq,DEBUG=True)
	#normalize(14,'GGCGGC','GGC','ATGATGGCGGCGGCGTAGTAG',0,DEBUG=True)

	#lvar1 = len(var1_ln.ref)
	#fullref = var1_ln.ref + refseq[lvar1+var1_ln.pos-offset:var1_rn.pos+lvar1-offset]
	#fullalt = var1_ln.alt + refseq[lvar1+var1_ln.pos-offset:var1_rn.pos+lvar1-offset]
	#print('var1-full',var1_ln.pos,fullref,fullalt)
