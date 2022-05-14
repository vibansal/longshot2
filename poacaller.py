from __future__ import print_function
import pysam,os,sys,glob,math,argparse,time,random
import subprocess
CIGAR_OP = ['M','I','D','N','S','H','P','=','X'];
import heapq
from pysam import VariantFile,TabixFile
import resource
from src.subread import SubRead
from src.errorbyread import read_error_stats
from src.interval import Interval, kmer_complexity
from src.printvariants import *
from src.bcftools import combinevcfs 

class VarCall:
	def __init__(self,args): #bamfile,fastafile,outdir,windowsize=150,vcf=None,filterreads=False):
		self.bam = args.bam
		self.readlist = []; ## this is a heap
		self.reads_heap=0;
		self.chrom=None; self.start=0; self.end=0; ## region in which to do variant calling

		self.params = args 

		self.windowsize=args.windowsize;
		self.flank = args.flank
		self.maxwindowsize = 3*args.windowsize
		self.minentropy = 2.0
		self.chromlength=0;
		self.variant_buffer = []
	
		self.open_files()

	def open_files(self):
		#if self.params.giab == None: self.pybed = None
		#else: self.pybed = TabixFile(self.params.giab) 

		if self.params.out!= None: 
			self.outfile=open(self.params.out,'w')
			self.logfile=open(self.params.out+'.log','w')
		else: 
			self.outfile=sys.stdout
			self.logfile=sys.stdout

		if self.params.vcf == None:	
			self.pyvcf = None
		else: 
			try:
				self.pyvcf = TabixFile(self.params.vcf)
			except OSError:
				print('index file not found for VCF, provide indexed vcf file',file=sys.stderr)
				#subprocess.call('bgzip ' +  self.params.vcf + ' > ' + self.params.vcf + '.gz',shell=True)
				#subprocess.call('tabix  ' +  self.params.vcf+'.gz',shell=True)
				#self.pyvcf = TabixFile(self.params.vcf+'.gz')
				self.pyvcf = None

		if self.outfile != sys.stdout: print_vcf_header(self.params.ref,self.outfile)
		if self.bam.endswith('cram') or self.bam.endswith('CRAM'): 
			self.bamreader = pysam.AlignmentFile(self.bam, "rc",reference_filename=self.params.ref)
		else: 
			self.bamreader = pysam.AlignmentFile(self.bam, "rb");
		self.pyFasta = pysam.Fastafile(self.params.ref);

	def close_filehandles(self):
		if self.outfile != sys.stdout: self.outfile.close()
		if self.pyvcf != None: self.pyvcf.close()
		#if self.pybed != None:self.pybed.close()
		self.pyFasta.close()
		self.bamreader.close();

	def parse_region(self,interval):
		self.chrom=interval.split(':')[0]; #print(self.chrom)
		self.chromlength = self.pyFasta.get_reference_length(self.chrom);
		if len(interval.split(':')) ==2: 
			if len(interval.split(':')[1].split('-')) ==2: 
				self.start = int(interval.split(':')[1].split('-')[0])
				self.end = int(interval.split(':')[1].split('-')[1])
				if self.end > self.chromlength: self.end = self.chromlength ## if it exceeds
			else: 
				p = int(interval.split(':')[1]);
				self.start = p-2000; self.end = p+2000;
		else: ## full chromosome
			self.start=0;
			self.end = self.chromlength;

	def process_bedfile(self,bedfile):
		with open(bedfile) as fh: regions = fh.readlines()
		for line in regions: 
			region = line.strip().split('\t')
			self.chrom = region[0]; self.start =int(region[1]); self.end = int(region[2])
			self.readlist = []
			self.reads_heap=0	
			self.process_bam(Sliding_Windows=False) ## entire region as one chunk
			

	def process_interval(self,interval,prev_right,INCLUDE_INPUT=True):
		s = max(0,interval.left-self.params.windowsize); e = (interval.right+self.params.windowsize);
		interval.paddedseq = self.pyFasta.fetch(self.chrom,s,e).upper();
		#interval.leftpad = self.pyFasta.fetch(self.chrom,s,interval.left).upper();
		#interval.rightpad = self.pyFasta.fetch(self.chrom,interval.right,e).upper();
		interval.paddedstart = s
		interval.populate_slices(self.readlist,interval.seq,outfile=self.logfile)
		interval.find_variants_POA(outfile=self.logfile);
		input_varlist = get_input_variants(self.pyvcf,self.chrom,prev_right,interval.right) ## extract longshot variants overlapping window from pyvcf TabixFile
		nvar = len(input_varlist)
		interval._print_(nvar,outfile=self.logfile)
		interval.combine_variants(outfile=self.outfile)
		self.variant_buffer +=interval.varlist 
		if INCLUDE_INPUT: self.variant_buffer +=input_varlist
		output_variants(self.variant_buffer,interval,outfile=self.outfile,logfile=self.logfile)
		print("#####################################################################################################################\n",file=self.logfile);

	def move_to_next_interval(self,left,right):
		left = right-self.flank ## ensures that start of window is also high-complexity
		right = left + self.params.windowsize;
		if right > self.end: right = self.end
		while right < self.end and right-left < self.maxwindowsize:
			ent,maxcount = kmer_complexity(self.pyFasta.fetch(self.chrom,right-self.flank,right),self.flank)
			if ent > self.minentropy and maxcount[0] < 5: break
			right += self.flank
			print('move boundary right',ent,right,file=self.logfile)
			if right-left >= self.maxwindowsize: print('max window size limit reached',file=self.logfile)
		return left,right


	## we can process two bams jointly as well, maintain two readlists
	def process_bam(self,DEBUG=False,Sliding_Windows=True):	
		self.reads_heap=0; processedreads=0;
		windows=0
		if Sliding_Windows: left = self.start; right = left + self.params.windowsize; 
		else:	left= self.start; right = self.end
		prev_right=left
		prev_left=left

		bam_iter= self.bamreader.fetch(self.chrom,self.start,self.end)
		is_bam_finished =False
		read_start = left
		while right < self.end+1 and not is_bam_finished: 
			##fetch new reads from bamfile until read.start exceeds interval.end/right
			while read_start < right: 
				try: read = next(bam_iter)
				except StopIteration: is_bam_finished = True; break
				if read.flag >= 256 or read.flag ==4 or read.mapping_quality <self.params.minMQ: pass
				else:
					subread = SubRead(read);
					subread.estats = read_error_stats(read,pflag=False,min_length=500);
					processedreads+=1;
					heapq.heappush(self.readlist,subread); 
					self.reads_heap +=1;
					read_start = read.reference_start

			if self.reads_heap > 0:
				interval = Interval(self.chrom,left,right,params=self.params,seq=self.pyFasta.fetch(self.chrom,left,right).upper())
				self.process_interval(interval,prev_right) ## variant calling

			windows +=1
			if windows%5000==0 and self.reads_heap >0: 
				print('processed',processedreads,'position',subread.refpos,file=sys.stderr)
				print('mem usage kb',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000,file=sys.stderr);

			prev_left = left
			prev_right =right
			if Sliding_Windows: left,right = self.move_to_next_interval(prev_left,prev_right)

			## delete those reads from heap for which read.reference_end < start of window
			while self.reads_heap > 0 and self.readlist[0].reference_end < left:
				if DEBUG: print('deleting read from heap',self.readlist[0].readid,self.readlist[0].reference_end,right,self.reads_heap); 
				heapq.heappop(self.readlist); 
				self.reads_heap -=1;
	
		output_variants(self.variant_buffer,interval,outfile=self.outfile,logfile=self.logfile,FINAL=True) ## any variants left in buffer
	
###############################################


def parseargs(empty=False):
    parser = argparse.ArgumentParser(description='## PROGRAM TO RUN SPOA on small windows of a BAM file, can be run on haplotype separared BAMs for variant discovery \n requirements: bcftools, bgzip and tabix \n \n')
    parser.add_argument('-w','--windowsize', nargs='?', default=300,type = int, help='window size for POA calculation')
    parser.add_argument('-m','--minreads', nargs='?', default=4,type = int, help='minumum reads required for POA consensus')
    parser.add_argument('--minMQ', nargs='?', default=10,type = int, help='minumum mapping quality of reads used for POA consensus')
    parser.add_argument('--maxdepth', nargs='?', default=500,type = int, help='maximum number of reads allowed for POA consensus')
    parser.add_argument('--flank', nargs='?', default=50,type = int, help='flanking sequence length for consensus, default 30 bp')
    parser.add_argument('--ref','--fasta', nargs='?', default=None,type = str, help='reference fasta file to which reads are mapped')
    parser.add_argument('-b','--bam', nargs='?', default=None,type = str, help='bam file')
    parser.add_argument('--no-hap',dest='nohaps',default=False,action='store_true',help='do not use haplotype tags from BAM file (if available) for POA')
    parser.add_argument('--addref',dest='addref',default=False,action='store_true',help='add reference sequence as first sequence in POA')
    parser.add_argument('--usecigar',dest='usecigar',default=True,action='store_true',help='add first read using cigar alignment in POA')
    parser.add_argument('--cluster',dest='cluster',default=True,action='store_true',help='group identical reads to speed up processing, for CCS reads')
    parser.add_argument('--compare',dest='compare',default=False,action='store_true',help='compare read counts supporting variant between two haplotypes')
    parser.add_argument('--filterreads', dest='filterreads',default=False,action='store_true',help='filter out sub-reads that match reference sequence')
    parser.add_argument('--pre',dest='prefilter',default=False,action='store_true',help='calculate per-read-error rates for POA ordering')
    parser.add_argument('-o','--out', nargs='?', default=None,type = str, help='output vcf file,combined')
    parser.add_argument('--vcf', nargs='?', default=None,type = str, help='VCF file with candidate variants, gzipped and indexed')
    #parser.add_argument('-V','--verbose', nargs='?', default=0,type = int, help='verbose print')
    parser.add_argument('-i','--region', nargs='?', default=None,type = str, help='target genomic interval')
    parser.add_argument('--bed', nargs='?', default=None,type = str, help='bedfile with intervals')

    if empty: return parser

    args = parser.parse_args()
    if len(sys.argv) < 4 or (args.bam == None) or args.ref == None or (args.region == None and args.bed == None): 
        parser.print_help()
        sys.exit(1)

    return args

######################################################################################################################


if __name__ == "__main__":
	args = parseargs(empty=False)
	varcaller = VarCall(args)
	if args.region != None: 
		print('options',args.bam,args.ref,args.region,args.out)
		varcaller.parse_region(args.region)
		varcaller.process_bam()
		#print('test mode... exiting',file=sys.stdout); sys.exit()
		varcaller.close_filehandles()
		if args.out != None: combinevcfs(args.out,args.ref,BCFTOOLS="/home/vbansal/Public/tools/bcftools-1.12/bcftools") 
	elif args.bed != None:
		varcaller.process_bedfile(args.bed)
		varcaller.close_filehandles()
	else: 
		print('specify a region or a bed file',file=sys.stderr)

