import os,sys
from .interval import Interval
from .variant import Variant,VCFrecord
from .variant_functions import *
from pysam import VariantFile,TabixFile,Fastafile

def print_vcf_header(fasta_file,outfile=sys.stdout):
	genome_fasta_open = Fastafile(fasta_file)
	#pyvcf = VariantFile(vcffile)
	#print(pyvcf.header,file=outfile,end='')
	#pyvcf.close()
	print("##fileformat=VCFv4.1",file=outfile)
	print("##source=POAvariantcaller",file=outfile)
	print("##reference="+fasta_file,file=outfile)
	print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",file=outfile)
	print("##FORMAT=<ID=RC,Number=1,Type=String,Description=\"coverage of haplotype\">",file=outfile)
	print("##FORMAT=<ID=RD,Number=1,Type=String,Description=\"total read depth\">",file=outfile)
	print("##FORMAT=<ID=WI,Number=1,Type=String,Description=\"Window for POA calling\">",file=outfile)
	print("##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Statistics for POA calling\">",file=outfile)
	print("##FORMAT=<ID=RT,Number=1,Type=Float,Description=\"fraction of reads supporting alt vs (ref+alt)\">",file=outfile)
	#print("##FORMAT=<ID=VIW,Number=1,Type=Float,Description=\"number of variants in window\">",file=outfile)
	print("##FORMAT=<ID=HP,Number=1,Type=String,Description=\"Haplotype for POA calling\">",file=outfile)
	#print("##FORMAT=<ID=SRC,Number=1,Type=String,Description=\"Source of the variant, input or POA-graph\">",file=outfile)
	print("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the variants, POA or input\">",file=outfile)
	print("##FILTER=<ID=dn,Description=\"In a dense cluster of variants\">",file=outfile)
	print("##FILTER=<ID=dp,Description=\"Exceeds maximum depth\">",file=outfile)
	print("##FILTER=<ID=sb,Description=\"Allelic strand bias\">",file=outfile)

	for chrom in genome_fasta_open.references:
		print("##contig=<ID="+chrom+",length="+str(genome_fasta_open.get_reference_length(chrom))+">",file=outfile)
	print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",file=outfile)
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
	genome_fasta_open.close()

	##contig=<ID=X,length=155270560>


def get_input_variants(pyvcf,chrom,start,end): ## store as a list of Variant objects (only vcf part is populated)
	if pyvcf==None: return []
	input_variants = []
	for rec in pyvcf.fetch(chrom,start,end):  ## tabix file
		v = rec.split('\t')
		if v[6] == 'RefCall': continue ## deepvariant filter 
		#print('inputv',v)
		varnew = Variant(pos=int(v[1]),ref=v[3],alt=v[4],source='caller');
		if len(varnew.ref) != len(varnew.alt): varnew.vtype = 'indel'
		else: varnew.vtype = 'snv'
		varnew.vcf = VCFrecord(chrom=v[0],pos=int(v[1]),identifier=v[2],ref=v[3],alt=v[4],qual=float(v[5]),filt=v[6],info=v[7].split(';'),formatfield=v[8],genotype=v[9])
		if ',' in v[4]: varnew.vcf.alts= v[4].split(',')
		else: varnew.vcf.alts = [varnew.vcf.alt]
		varnew.Wstart = start
		varnew.Wend = end
		varnew.vcf.leftnormal= Vartuple(pos=varnew.vcf.pos,ref=varnew.vcf.ref,alt=varnew.vcf.alt,alleles=[])
		input_variants.append(varnew) 
	return input_variants

def output_variants(variant_buffer,interval,outfile=sys.stdout,logfile=sys.stdout,DEBUG=True,FINAL=False,normalize=False):
	variant_buffer.sort(reverse=True)
	n = len(variant_buffer)
	while n > 0:
		var = variant_buffer.pop()
		if var.vcf.pos < interval.left or FINAL: 
			if not var.is_duplicate:# and var.flag != 'filt':  ## changed this for CCS
				var.print_VCFformat(outfile=outfile,logfile=logfile)
			n -=1
		else:
			variant_buffer.append(var)
			break

