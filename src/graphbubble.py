import sys
import statistics
from .variant import Variant

"""
functions to extract bubbles(variants) from two paths (consensus,reference) derived from POA graph 
each path has a sequence and edge weights
"""

class CPath: ## compressed sequence path, can also represents a bubble (two paths that start and end at same node in POA graph)
	def __init__(self,position,is_bubble=False):
		self.is_bubble=is_bubble
		self.position = position
		self.length=0
		self.sequence = []
		self.counts = []
		self.vartype= None
		self.used=False ## consumed during bubble creation

		#self.var = None ## variant record

		self.dprev=0; self.dnext=0 ## distance to previous and next bubbles
		self.prevbase = ''
		
		self.length_r=0
		self.sequence_r = []
		self.counts_r = []

		mincov = 100; mincov_ref = 100; fraction = 0.0; ratio_alleles=0.0

	def calc_stats(self,reads):
		self.mincov = min(self.counts)
		self.mincov_ref = min(self.counts_r)
		self.ratio_alleles = float(self.mincov)/((self.mincov+self.mincov_ref))
		self.fraction = float(self.mincov)/reads;

	def _print_(self,outfile=sys.stdout):
		if self.is_bubble: 
			self.print_bubble(outfile=outfile)
		else:
			mincov = min(self.counts)
			maxcov = max(self.counts)
			l = self.length
			if self.length >= 20: print('path',self.length,''.join(self.sequence[0:min(l,20)] + ['.','.'] + self.sequence[max(20,l-20):l]),'cov',mincov,maxcov,file=outfile)
			else: print('path',self.length,''.join(self.sequence),','.join([str(c) for c in self.counts]),file=outfile)

	def print_bubble(self,outfile=sys.stdout):
		print(self.dprev,self.dnext,'---------bubble...',self.position+1,''.join(self.sequence_r),''.join(self.sequence),end=' ',file=outfile)
		print(self.mincov_ref,self.mincov,'stats',round(self.fraction,3),round(self.ratio_alleles,3),end='\t',file=outfile)
		print('cov-ref',','.join([str(c) for c in self.counts_r]),'alt',','.join([str(c) for c in self.counts]),file=outfile)
		#print(self.var.pos,self.var.ref,self.var.alt,self.var.flag,self.var.vtype,'hap',self.var.haptag,end=' ')

	def bubble_to_variant(self,chrom=None,haptag='0',maxcov=0): ## convert CPath bubble object to a 'Variant' object 
		lref = len(self.sequence_r)
		lalt = len(self.sequence)
		seq_ref = ''.join(self.sequence_r) ## trim one base off end, shared for both alleles
		seq_alt = ''.join(self.sequence)
		if seq_ref == seq_alt: flag = 'equal'
		delta = lalt-lref
		if delta == 0: 
			vtype ='SNV'
			POS = self.position +1
			REF = seq_ref; ALT = seq_alt
		else: 
			vtype = 'INDEL'+':'+str(delta)
			POS = self.position + 1
			#REF = self.prevbase + seq_ref; ALT = self.prevbase + seq_alt
			REF = seq_ref; ALT = seq_alt
		var= Variant(chrom=chrom,vtype=vtype,pos=POS,ref=REF,alt=ALT,source='POA',haptag=haptag,maxcov=maxcov)
		var.mincov = self.mincov
		var.mincov_ref = self.mincov_ref
		var.ratio_alleles = self.ratio_alleles
		var.fraction = self.fraction
		return var

	def merge_bubbles(self,other,middle):
		self.sequence += middle.sequence + other.sequence
		self.counts += middle.counts + other.counts
		self.counts_r += middle.counts + other.counts_r
		self.sequence_r += middle.sequence + other.sequence_r
		self.length_r += middle.length + other.length_r
		self.length += middle.length + other.length

	def check_bubbles(self,bubble2,middle,DEBUG=True,outfile=sys.stdout):
		ref = self.sequence_r +middle.sequence + bubble2.sequence_r
		alt = self.sequence + middle.sequence + bubble2.sequence
		i = len(ref); j = len(alt)
		while i > 0 and j >0 and ref[i-1] == alt[j-1]:
			i -=1; j -=1
		s=0
		while s < i and s < j and ref[s] == alt[s]: s +=1
		if DEBUG: print(''.join(ref),''.join(alt),'-->',''.join(ref[s:i])+ '|' + ''.join(alt[s:j]),end=' ',file=outfile)

		if self.vartype == 'SNV' or bubble2.vartype == 'SNV': ## if one of the two bubbles is a SNV, only merge if the combined variant is a clean indel
			if (i-s ==0 or j-s==0):
				if DEBUG: print('tomerge-simple',file=outfile);
				return 1	
			else: 
				if DEBUG: print('nonmerge',file=outfile)
				return 0

		if i-s == j-s and i-s <=2: 
			if DEBUG: print('tomerge-snv',file=outfile)
			return 1
		elif i-s <2 or j-s <2: 
			if DEBUG: print('tomerge-indel',file=outfile)
			return 1
		else: 
			if DEBUG: print('nonmerge',file=outfile)
			return 0


#####################################################################################################################

class Subgraph:
	def __init__(self,cons_path,ref_path,outfile=sys.stdout):
		self.nodes = {} 	## nodes is a sub-graph represented using a dictionary, key is 'node-id'
		self.is_ref = {}
		self.first = -1
		self.last = -1
		self.diff=0
		self.trimmed=0
		self.pathlist = [] ## array of CPath objects 
		self.nbubbles=0
		self.refbases_covered=0 
		self.outfile=outfile
		self.add_path(cons_path,label='alt') ## add alt first
		self.add_path(ref_path,label='ref')
		self.min_ends = 6 ## first(last) path in graph should be of length >= this value
		self.bubble_list = []

	def add_path(self,graph_path,label='ref',mincov=1): ## graph_path = list of edges where edge = (s,e,wt,base-label-incoming,position_for_ref)
		ln = len(graph_path)
		startp = 0
		endp = ln
		if label == 'alt': ## trim '0' weight edges from beginning and end of path due to reference dummy
			while startp < ln and graph_path[startp][2] ==0: startp +=1
			while endp > startp and graph_path[endp-1][2] ==0: endp -=1
			if startp > 0 or endp < ln-1: 
				print('removed dummy ref edges from path',startp,endp,ln,file=self.outfile)
				self.trimmed=1

		for index in range(startp,endp):
			edge = graph_path[index]
			start =edge[0]; end= edge[1]; wt = edge[2]
			if label == 'ref': self.is_ref[start] = True; self.is_ref[end] = True
			if self.first<0 and edge[2] >= mincov and label == 'ref': self.first = start; 
			if edge[2] >= 1 and label == 'ref': self.last = end
			## add to graph
			try:
				#if self.nodes[start][0][0] == end: ## same edge 
				self.nodes[start].append((edge[1],edge[3],edge[2],label,edge[4]))
				if self.nodes[start][0][0] != self.nodes[start][1][0]: self.diff +=1
			except KeyError: 
				self.nodes[start] = [(edge[1],edge[3],edge[2],label,edge[4])]; 
				if label == 'ref': self.diff +=1

	def extract_bubbles_graph(self,ref_start,reads): ## main function
		if self.diff == 0: return [],0
		else:
			self.find_bubbles(ref_start)
			self.bubble_list= self.clean_bubbles(reads)
			self.nbubbles =len(self.bubble_list)
			print('count',self.nbubbles,'refbases',self.refbases_covered,'\n',file=self.outfile)

	"""
	go over the list of paths, keep adding to new list, if we encounter closeby, merge and add, continue... only one iteration is sufficient
	separation could be more than 20 for low-complexity indels...
	"""
	def clean_bubbles(self,reads,merge_distance=20,MINCOV=3,DEBUG=True):
		prev_bubble = -1
		dist_last_bubble=1000
		newpathlist = []
		first_path=-1; last_path=-1
		np=len(self.pathlist)
		self.refbases_covered=0

		## trim paths from beginning and end until we find a shared path of length >= min_ends
		for i in range(np):
			path = self.pathlist[i]
			if not path.is_bubble and first_path < 0 and path.length >= self.min_ends: first_path = i
			if not path.is_bubble and path.length >= self.min_ends: last_path = i+1

		if first_path > 0 or last_path < np: print('trimmed paths',first_path,last_path,np,file=self.outfile)

		"""
			logical flaw: self.pathlist[i-1] used as middle sequence incorrectly....when we have Bubble1-middle-Bubble2-Bubble3 
			fixed 02/05/22
		"""

		dummy_path = CPath(0) ## empty path for merging
		for i in range(first_path,last_path): self.pathlist[i].used = False
		for i in range(first_path,last_path):
			path = self.pathlist[i]
			if not path.is_bubble:
				self.refbases_covered += path.length
				if i > first_path: self.pathlist[i-1].dnext = path.length
				if DEBUG: path._print_(outfile=self.outfile)
				
			else:
				path.dprev = self.pathlist[i-1].length
				path.prevbase = self.pathlist[i-1].sequence[-1] ## last base of previous sequence
				path.calc_stats(reads)
				self.refbases_covered += path.length_r
				flag=0
				path_mincov = min(path.counts)
				if DEBUG: path._print_(outfile=self.outfile)
				if prev_bubble >=0: dist_last_bubble = self.pathlist[i].position-self.pathlist[prev_bubble].position-self.pathlist[prev_bubble].length

				## check if this bubble and prev bubble need to be merged
				if prev_bubble >=0 and dist_last_bubble  < merge_distance and (self.pathlist[prev_bubble].vartype == 'INDEL' or path.vartype=='INDEL'):
					if path_mincov >=MINCOV:
						if prev_bubble == i-1: flag =self.pathlist[prev_bubble].check_bubbles(self.pathlist[i],dummy_path,outfile=self.outfile) 
						else: flag = self.pathlist[prev_bubble].check_bubbles(self.pathlist[i],self.pathlist[i-1],outfile=self.outfile);

				if flag ==1: ## merge_bubbles, there could be small path between bubbles, 'middle-path=self.pathlist[i-1]'
					if prev_bubble == i-1 or self.pathlist[i-1].used == True: 
						flag =self.pathlist[prev_bubble].merge_bubbles(self.pathlist[i],dummy_path) 
						self.pathlist[i].used = True
					else: 
						flag = self.pathlist[prev_bubble].merge_bubbles(self.pathlist[i],self.pathlist[i-1]);
						self.pathlist[i].used = True; self.pathlist[i-1].used = True
					self.pathlist[prev_bubble].calc_stats(reads) ## updated
				else: 
					newpathlist.append(path)
					prev_bubble = i
				#print('debug',i,prev_bubble,dist_last_bubble)

		return newpathlist

	## takes the graph dictionary from merge_paths and generate pathlist (sequence of simple-paths and bubbles)
	def find_bubbles(self,ref_offset):
		pathlist = []
		flag =1
		s =self.first
		while s in self.nodes:
			edges = self.nodes[s]
			num_edges = len(edges)
			#print('debug bubble finding',s,edges)
			if num_edges ==1: ## new case where this path is only in reference
				s = edges[0][0]
			elif num_edges ==2 and edges[0][0] == edges[1][0]: ## shared edge in both paths
				if flag==1:
					#print('debug start',ref_offset,edges[1][4]+start) 
					pathlist.append(CPath(ref_offset+edges[1][4],is_bubble=False))
					flag=0
				pathlist[-1].length +=1
				pathlist[-1].sequence.append(edges[0][1])
				pathlist[-1].counts.append(edges[0][2])

				s =edges[0][0]

			elif num_edges ==2: ## divergence 
				s1,edgelist1,s2,edgelist2,terminate_flag = self.get_bubble_edges(s)
				#print('bubble',s,'|',s1,edgelist1,s2,'ref',edgelist2)
				l2 = len(edgelist2)
				if terminate_flag and l2>=10: print('terminate flag',l2,s1,s2,file=self.outfile)
				if s1 != s2: # unable to align ref/alt paths 

					break
				else:
					pathlist.append(CPath(ref_offset+edgelist2[0][4],is_bubble=True))
					flag = 1
					for edge in edgelist1: pathlist[-1].sequence.append(edge[1]); pathlist[-1].counts.append(edge[2]); pathlist[-1].length +=1
					for edge in edgelist2: 
						pathlist[-1].sequence_r.append(edge[1]); pathlist[-1].counts_r.append(edge[2]); 
						pathlist[-1].length_r +=1
					if pathlist[-1].length != pathlist[-1].length_r: pathlist[-1].vartype = 'INDEL'
					else: pathlist[-1].vartype = 'SNV'
				s = s2
				if s not in self.is_ref: break
			else: break
		self.pathlist = pathlist


	def get_bubble_edges(self,s):
		## start from 's' and find path of non-ref nodes in graph that ends at a 'ref' node 'e'
		edgelist1 = []
		curr_edge = self.nodes[s]
		terminate_flag=False
		next_node = -1
		while True:
			for neighbor in curr_edge: 
				if neighbor[3] == 'alt': edgelist1.append(neighbor); next_node = neighbor[0]

			if next_node in self.is_ref: break
			try: curr_edge = self.nodes[next_node]
			except KeyError: 
				terminate_flag=True
				break ## not in graph
		last_node = next_node

		## find the path of 'ref-nodes' from s to 'e'
		edgelist2 = []
		curr_edge = self.nodes[s]
		while True:
			for neighbor in curr_edge: 
				if neighbor[3] == 'ref': edgelist2.append(neighbor); next_node = neighbor[0]
			if next_node == last_node and not terminate_flag: break
			try: curr_edge = self.nodes[next_node]
			except KeyError: break

		return last_node,edgelist1,next_node,edgelist2,terminate_flag

	"""
	## add reference path that is not supported by consensus, not a good idea since it is not aligned...
	if len(edgelist2) >= 2: 
		pathlist.append(CPath(ref_offset,is_bubble=False))
		pathlist[-1].length =len(edgelist2)
		for edge in edgelist2:
			pathlist[-1].sequence.append(edge[1])
			pathlist[-1].counts.append(edge[2])
		print('adding terminal ref-path',s1,s2,len(edgelist2)); 
	"""
