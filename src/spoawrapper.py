import math,argparse,time,random,os,sys
#sys.path.append(os.getcwd()+'/spoa_python/')
#from .spoa_python import spoa
import spoapy
from collections import Counter
#from graphbubble import merge_paths,find_bubbles

"""
code to take a set of sequences and run SPOA
for CCS reads, heaviest sequence should be added first, sort by sequence... for CLR, add in same order
add consensus from 2nd haplotype to graph, same alignment frame, easier to find homozygous variants vs hets... 

single graph with hap1-cons, hap1-ref, hap2-cons, hap2-ref (paths, bubbles, counts) 

"""
	
def analyze_msa(graph):
	msab = graph.generate_msa();
	msa = []
	for seq in msab: msa.append(seq.decode('UTF-8'))
	nseq = len(msa)
	cols = len(msa[0]) 
	total =0.0
	deltas= []
	for seq in msab: deltas.append(0.0);
	for col in range(cols):
		counts = {'A':0,'C':0,'G':0,'T':0,'-':0}
		delta = {'A':0.0,'C':0.0,'G':0.0,'T':0.0,'-':0.0}
		entropy = 0.0
		for row in range(nseq): counts[msa[row][col]] +=1
		for c in counts.keys(): 
			if counts[c] ==0: continue
			entropy += float(counts[c])*math.log(counts[c])
		total += entropy
		for c in counts.keys(): 
			if counts[c] ==0: original = 0.0
			else: original = float(counts[c])*math.log(counts[c])
			if counts[c] < 2: new1 = 0.0;
			else: new1= float(counts[c]-1)*math.log(counts[c]-1)
			delta[c] = new1-original
		for row in range(nseq):
			deltas[row] += delta[msa[row][col]]
		"""
		"""
	for row in range(nseq): print(round(deltas[row],2),msa[row])
	print('entropy',cols,total)

def query_sequence_graph(graph,alignment,sequence,delta=0,output='tuple',DEBUG=False):
	l = len(alignment)

	## need this when using OV mode for graph alignmnet
	#if alignment[l-1][1] < 0: print('ERROR',file=sys.stdout)
	#while alignment[l-1][1] < 0: l -=1 ## temporary fix... 02/07/22
	nodelist = [alignment[i][0] for i in range(l)] ## graph nodes in path
	weightlist = graph.get_path_weights(l,nodelist)
	#if output == 'vector': 
	cov= [int(weightlist[i]/2) for i in range(l-1)];
	if DEBUG: print(sequence,len(alignment),len(sequence),'\n',alignment)
	path = [[alignment[i][0],alignment[i+1][0],int(weightlist[i]/2),sequence[i+1],alignment[i+1][1]] for i in range(l-1)]
	#print('ref',cov)
	return cov,path


## scale edge weights of reference path by 0.5 before calling consensus_edges(), will prefer non-ref path but not bogus ones...
def reduce_reference_consensus(graph,alignment,refcov,scalefactor=0.5):
	prior_weights = {}
	j=0
	for i in range(len(alignment)-1):
		sn = alignment[i][0]; en = alignment[i+1][0]; 
		if sn >= 0 and en >=0: 
			prior_weights[(sn,en)] = refcov[j]
			if refcov[j] >= 2: wt= graph.update_edge(sn,en,int(-0.5*refcov[j])); 
		j +=1

	svec,Edges = graph.consensus_edges();
	cons_seq = svec.decode('UTF-8')
	cvec = []
	for edge in Edges: ## if edge matches ref-edge, need to get correct weight 
		try: 
			w= prior_weights[(edge.start_node_id(),edge.end_node_id())]
			cvec.append(w)
		except KeyError: cvec.append(int(edge.weight()/2))
	cvec.append(0)
	return cons_seq,cvec,Edges

def add_cigar_graph(graph,refseq,ref_alignment,sequence,cigarlist,outfile=sys.stdout): ##
	#assumption is that refseq and refseq cigar is correct, check it before adding, can be done during loop
	#print(refseq,ref_alignment,sequence,cigarlist)
	"""
	alignment of sequence [1...n] to graph generates list of two tuples (graph-node,sequence-index) 
	convert cigar to graph alignment object: 
	1. deletion -> simple, (j,i) followed by (j+20,i+1)  
	2. mismatch -> (j,i) (j+1,i+1), (j+2,i+2)  | same as match | don't need extended cigar 
	3. insertion -> (j,i)... (-1,i+1) (-1,i+2) ..... (j+1,i+l)

	alignment object can be constructed from list of 2-tuples: AL = spoapy.Alignment.from_pairs([(1,1),(2,2)])
	once alignment object is there, do graph.add_alignment(AL,sequence,1) 
	"""
	cl = len(cigarlist)
	if cl < 2: 
		print('no need',cigarlist,file=outfile)
		return 0

	pairs = []
	ref_index =0; seq_index=0; mismatches=0
	FLAG=0
	for i in range(cl):
		cigar = cigarlist[i]
		l = len(cigar)
		length = int(cigar[0:l-1])
		if length <0:
			FLAG=1
		op = cigar[l-1]
		if op == 'M' or op == '=' or op == 'X': 
			for j in range(length): 
				pairs.append((ref_index,seq_index)); ref_index+=1; seq_index +=1
		elif op == 'D': 
			ref_index += length
			if i ==0 or i==cl-1: FLAG=1
		elif op == 'I': 
			for j in range(length): pairs.append((-1,seq_index)); seq_index +=1
			if i ==0 or i==cl-1: FLAG=1

	if FLAG ==0 and seq_index == len(sequence) and ref_index == len(refseq):
		AL = spoapy.Alignment.from_pairs(pairs)
		#print(pairs,'alignment')
		graph.add_alignment(AL,sequence,0)
		print('adding first read to graph via cigar...',cigarlist,file=outfile)
		return 1
	else:
		print(cigarlist,'ERROR not adding',FLAG,seq_index,len(sequence),ref_index,len(refseq),file=outfile)
	return 0	
	#print('cigarpairs',pairs)
	#print(refseq,ref_alignment,sep='\n')
	#sys.exit()

	"""	
	newseq = refseq[0:5] + 'C' + refseq[6:50]; print(refseq[0:50]); print(newseq[0:50])
	"""


## call SPOA library (C++ with python wrapper) to get consensus and weights vector cvec
def spoa_python(sequences,match=5,mismatch=-4,gapO=-8,gapE=-6,mincov=3,weighted=True,refseq=None,ADD_FIRST=True,DIPLOID=False,align_type='SW',CIGAR_FIRST=True,outfile=sys.stdout):
#def spoa_python(sequences,match=1,mismatch=-4,gapO=-6,gapE=-5,mincov=3,weighted=True,refseq=None,ADD_FIRST=True,DIPLOID=False,align_type='SW',CIGAR_FIRST=True):


	if align_type == 'SW':
		engine = spoapy.AlignmentEngine(spoapy.AlignmentType.SW,m=match,n=mismatch,g=gapO,e=gapE) ## can add options m=5, n=-4, g=-8, e=-6, q=-10, c= -4
	elif align_type =='NW':
		engine = spoapy.AlignmentEngine(spoapy.AlignmentType.NW,m=match,n=mismatch,g=gapO,e=gapE) ## can add options m=5, n=-4, g=-8, e=-6, q=-10, c= -4
	elif align_type =='OV':
		engine = spoapy.AlignmentEngine(spoapy.AlignmentType.OV,m=match,n=mismatch,g=gapO,e=gapE) ## can add options m=5, n=-4, g=-8, e=-6, q=-10, c= -4
		
	graph = spoapy.Graph()
	first_added=0
	if refseq!= None and (CIGAR_FIRST or ADD_FIRST): 
		graph.add_sequence(engine,refseq,0); ## add refseq first
		ref_alignment = list(graph.align(engine,refseq)); ## returns two-tuple that can have -1 values for both...
		if CIGAR_FIRST: 
			first_added = add_cigar_graph(graph,refseq,ref_alignment,sequences[0][1],sequences[0][-1],outfile=outfile)
			#print(list(graph.align(engine,sequences[0][1])),sequences[0][1]); ## returns two-tuple that can have -1 values for both...

	if weighted:
		#ops =0; total=0 
		i=0
		for read in sequences:  
			#wt = read[0]
			#if first_added ==1 and i==0: wt = read[0]-1
			graph.add_sequence(engine,read[1],read[0]); 
			#i +=1
			#ops +=1; total+= read[1]
		#print('saved POA',total-ops,total)
	else:
		i=0
		for read in sequences: 
			#if first_added==1 and i==0: pass
			graph.add_sequence(engine,read[1],1) ## use order provided
			#i +=1

	if not ADD_FIRST and refseq != None and not CIGAR_FIRST:
		## create separate alignment engine for reference seq (NW rather than SW) 
		#engine_new = spoapy.AlignmentEngine(spoapy.AlignmentType.OV,m=match,n=mismatch,g=gapO,e=gapE) ## can add options m=5, n=-4, g=-8, e=-6, q=-10, c= -4
		graph.add_sequence(engine,refseq,0); ## add refseq last
		ref_alignment = list(graph.align(engine,refseq)); ## returns two-tuple that can have -1 values for both...

	if DIPLOID == False or refseq==None:
		svec,Edges = graph.consensus_edges();
		cons_seq = svec.decode('UTF-8')
		cvec =  [int(edge.weight()/2) for edge in Edges]; cvec.append(0)
		cons_path = [[Edges[i].start_node_id(),Edges[i].end_node_id(),int(Edges[i].weight()/2),cons_seq[i+1],-1] for i in range(len(Edges))]

	#analyze_msa(graph) 
	#graph.print_dot('vikas'); print('outputting graph for poa',refseq,file=sys.stderr); sys.exit()

	if refseq != None:
		#print('debug',len(ref_alignment),len(refseq),ref_alignment)
		ref_cov,ref_path = query_sequence_graph(graph,ref_alignment,refseq,output='vector'); ref_cov.append(0)
		#print(ref_alignment,'\n',list(graph.align(engine,cons_seq)))
		#print(refseq,'\n',ref_cov,'\n',cons_seq,'\n',cvec)
		#print('ref',len(ref_path),'\nalt',len(cons_path))
		
	else: ref_cov = None 

	if DIPLOID == True: 
		cons_seq,cvec,Edges = reduce_reference_consensus(graph,ref_alignment,ref_cov)
		cons_path = [[Edges[i].start_node_id(),Edges[i].end_node_id(),cvec[i],cons_seq[i+1],-1] for i in range(len(Edges))]

	#nodes,first,offset = merge_paths(cons_path,ref_path)
	#if first >=0: find_bubbles(nodes,first,offset,max(cvec))

	return cons_seq,cvec,ref_cov,(cons_path,ref_path)


#graph.update_edge(Edges[0].start_node_id(),Edges[0].end_node_id(),-2);
