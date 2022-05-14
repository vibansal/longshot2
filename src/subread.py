from __future__ import print_function
import pysam,os,sys
CIGAR_OP = ['M','I','D','N','S','H','P','=','X'];

class SubRead: ## stores bam read information relevant for slicing from full aligned read 
	def __init__(self,bamread):
		self.readid = bamread.query_name;
		self.cigartuples = bamread.cigartuples;
		self.ncigs = len(self.cigartuples);
		self.readseq = bamread.query_sequence;
		self.readlen = len(self.readseq)
		self.flag =bamread.flag;
		self.MQ = bamread.mapping_quality
		self.partial= False
		self.strand = '+'
		if bamread.flag & 16 == 16: self.strand = '-'
		try: 
			self.haptag = bamread.get_tag('HP'); ## haplotype tag for longshot or any tool, 0 is unassigned, 1/2 are assigned
			#if self.haptag ==0: self.haptag = -1;
		except KeyError: self.haptag = 0
		#BAM file with HP:i:1 or HP:i:2 depending on which haplotype, and also adds a PS tag for haplotyp block
		self.estats = [0.01,0.01,0.01,100]; ## error rate statistics useful for POA
		self.subcigar=None;
		self.subread = None;
		self.indelcount = [0,0] #(deletion,insertion)
		self.indexes = [-1,-1];
		self.refpos = bamread.reference_start;  ## updated as we move through the read
		self.readpos = 0;
		self.reference_start = bamread.reference_start;
		self.reference_end = bamread.reference_end;
		self.current = 0; # index in cigar that points to previous starting point, use to increase speed
		self.last=0;

	## for short reads, read could end in soft clip that should be included in window, not true for long reads

	## start and end of sub-read have to overlap with 'M/=/X/D' cigar operation
	## currently, function will only output subread when read overlaps entire window (start,end) | this doesn't work for short reads
	def getslice(self,start,end,DEBUG=False,min_overlap=-1,outfile=sys.stdout): ## bamread, self.current is a variable that is updated in each call to function
		cigar = []; # for that window
		indexes = [-1,-1]; ## start and end of slice on readsequence
		self.indelcount[0] = self.indelcount[1] = 0
		self.partial=False
		if self.reference_start > start or self.reference_end < end: ## if read starts after the window or ends, partial overlap
			self.partial = True
			#print('partial read',start,end,'read-s-e',self.refpos,self.reference_start,self.reference_end,self.readid,self.subcigar,indexes)
			self.subread = None; 
			if min_overlap <=0:return 0;
		refpos = self.refpos; readpos = self.readpos;
		flag =0;
		i = self.current; ## current cigar index
		while i < self.ncigs:
			ct = self.cigartuples[i];

			## move self.current forward to start of this window 
			if refpos+ct[1] > start and ct[0] != 4 and flag ==0: 
				self.refpos = refpos; self.readpos = readpos;
				self.current=i; flag=1

			if ct[0] == 0 or ct[0] == 7 or ct[0] ==8:  ## M/X/=
				if (refpos+ct[1] > start and refpos < end): ## aligned bases overlap interval
					off1 = 0; off2=0;
					if refpos <= start: off1 = start-refpos; indexes[0] = readpos + off1; ## only if refpos is less than start 
					if self.reference_start > start and indexes[0] < 0: indexes[0] = readpos ## partial
					if refpos + ct[1] >= end: 
						off2 = refpos+ct[1]-end; indexes[1] = readpos + end-refpos; ## only if read covers end
					if self.reference_end < end: indexes[1] = readpos
					lop = ct[1]-off1-off2;
					cigar.append(str(lop)+CIGAR_OP[ct[0]]);
					#print('stats',refpos,ct[1],start,end,lop,cigar[-1]);
				refpos += ct[1]; readpos += ct[1]; 
				if refpos > end + 5: 
					self.last=i;
					break;
			elif ct[0] == 1: ## insertion
				if refpos >= start and refpos < end: cigar.append(str(ct[1])+CIGAR_OP[ct[0]]); 
				self.indelcount[1] +=1
				readpos += ct[1]; 
			elif ct[0] == 2: ## deletion in reference, no base in read used
				if (refpos+ct[1] >= start and refpos < end): 
					off1 = 0; off2=0;
					if refpos <= start: off1 = start-refpos; indexes[0] = readpos
					elif refpos + ct[1] >= end: ## deletion crosses boundary of window 
						#off2 = refpos+ct[1]-end; 
						cigar.append(str(ct[1]+refpos-end-1)+CIGAR_OP[2])
						indexes[1] = readpos
						refpos += ct[1];  
						self.last=i+1
						break;
					else:
						cigar.append(str(ct[1]-off1)+CIGAR_OP[2]);
						self.indelcount[0] +=1
				refpos += ct[1];  
			elif ct[0] == 4: ## soft clip 
				readpos += ct[1]; 
			#if refpos > end: break; ## stop 
			i +=1;
		#self.last=i;
		self.subcigar=cigar;

		if self.partial and indexes[0] >=0 and indexes[1] >=0 and indexes[1]-indexes[0] >= min_overlap:
			self.subread = self.readseq[indexes[0]:indexes[1]];
			print('partial',self.readid,start,end,indexes,len(self.readseq),'seq',self.reference_start,self.reference_end,readpos)

		if indexes[1]-indexes[0] > 1 and indexes[0] >= 0:
			if indexes[1] < self.readlen+1: self.subread = self.readseq[indexes[0]:indexes[1]];
			else: 
				print('BUG',indexes[1],'exceeds read length',self.readlen,self.readid,self.cigartuples,self.subcigar,file=outfile)
				self.subread = None
		else: self.subread =None;

		self.indexes = [indexes[0],indexes[1]];
		if DEBUG: print('read',self.readid,start,end,indexes,len(self.readseq),'seq',self.subread,self.subcigar,self.current)
		return 1;

	# functions for heap comparison
	def __eq__(self,other):
		return self.reference_end == other.reference_end
	def __lt__(self,other):
		return self.reference_end < other.reference_end

	def printslice(self):
		pass;


