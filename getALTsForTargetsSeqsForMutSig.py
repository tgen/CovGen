#!/usr/bin/env python

# Copyright (C) 2017 Translational Genomics Research Institute

# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

# Major Contributors: Christophe Legendre

"""AIM: Create Three VCF files with Alternates mutations to each target sequences in provided bed File ; 
These VCF files  will further be annotated by snpEff (or any vcf annotator) to get the Effect and 
use that effect to create the coverage file for MutSig \nWARNING_REQUIREMENTS: The input BED file CANNOT have 
overlapping loci; If this is the case, please merge these loci."""


# Note: for a 3GB reference genome and a 3.5MB bedfile, the Virtual Memory usage will go up to 14GB of RAM
# the bottleneck of that script is the I/O step at the end when writing the 3 vcfs

try:
	import logging
	from multiprocessing import Process
	import multiprocessing
	from functools import partial
	from Bio import SeqIO
	from Bio.Seq import Seq
	from os.path import basename, exists
	from sys import exit, argv, stdout
	from natsort import natsorted, ns
	from datetime import date
except ImportError as err:
	print(err)
	raise ImportError('Import Error; Aborting')

__author__ = 'Christophe Legendre'
__company__ = "TGen"
__email__ = "clegendre@tgen.org"
__version__ = '0.1'
__copyright__ = """MIT License
	
	Copyright (c) 2017 Translational Genomics Research Institute
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE."""

#definition_functions
def file_check(parser, arg):
	if not exists(arg):
		parser.error("The file {0} does not exist!\nCheck your input".format(arg))
	else:
		return str(arg)

#returns the adequate length of the nucleotide of the 'input' arguments
def split_by_size(input, length):
	''' This is a very tricky function; it returns the trinucletide in the iput sequence, and take care of the
	 end of the sequence if the remaining is not the lenght of the expected tri-nucleotide'''
	return [input[i:i + length] for i in range(0, len(input), 1) if len(input[i:i + length]) == length]

def get_possible_alts(x):
	''' switch/case statement equivalent in Python'''
	# We need to manage the N or any letter that is not A, C, G or T  # TODO
	return {
		'A': ("C","G","T"),
		'C': ("A","G","T"),
		'G': ("A","C","T"),
		'T': ("A","C","G"),
		'N': ("N")
	}[x]

def get_ALTs(input):
	""" this time the lenght is one for the triNucleotide ; Why one ? because we want ALT mutation
		and therefore, we need to return the three others bases than the one at that locations
		We need to Check if we have N's as well and return None.
		input is the target sequence
		corrd is the contig name and starting position of the input target sequence such as: Y:123456
		lenght is the length of the sub-sequence, here it should be 1
	"""
	length=1 # this time the sliding window is not three as in the function split_by_size, but 1 as we want to go over 1 base at a time
	return( [ get_possible_alts(input[i:i+length]) for i in range(1, len(input)-1, length) if len(input[i:i + length]) == length ] )


def process_seq_for_alts(seq, seq_coords):
	""" the seq is actually one of the target sequence captured form the BED file ; 
		Here, we captured the trinucleotides for each frame  
	"""
	if seq_coords == "": return (None)
	list_sub_seqs = []
	if seq is None or seq == "": return (None)
	## ''' So here we are getting the information for one vcf per reading frame ; to do so we need information about the chr and pos as well '''
	for frame in range(1):
		## the number of Frame is one here as we slide one base at a time (not a trinucleotide)
		## this is the split_by_size funciton that does the sliding windows; this is the trick.
		seq_items = [ seq_coords, [ str(y) for y in split_by_size(seq[frame:], 3)],[x for x in get_ALTs(seq[frame:])] ]
		list_sub_seqs = list_sub_seqs + seq_items
	return (list_sub_seqs)
# note the list returned here is a list of lists of 3 lists and .. the items in the 3 lists are linked or aka they map to each other
# example: [ list1 is like [ 10:123456-123499 ], list2 is like [ATG, TGA, GAC,  ...], list3 is like [(A,C,G),(A,C,T),(C,G,T), ...]
# number of codons should be equal to the number of tuples in list3, and the trick is that the numbr of bases covered by the loci in list1 should
# be equal to the number of codons in L2 or the number of tuples in L3

def mergeLoci(bedf):
	import pybedtools
	logger.info("Merging Loci in BedFile ...")
	b = pybedtools.BedTool(bedf)
	b_merged = b.sort().merge()
	newBedFilename = bedf+".merged.vcf"
	b_merged.moveto(newBedFilename)
	logger.info("Merging Loci in BedFile ... DONE")
	return newBedFilename

def write2vcf(fh, pool_list_seqs, ALT_BASE_NUMBER):
	for list_seqs in pool_list_seqs:
		for i in range(len(list_seqs[1])):
			bases_in_codon = list(str(list_seqs[1][i]).upper())  ## we get the codon all in capital letters
			if "N" in bases_in_codon: continue
			## column ID could contain the Codon such as : str(list_seqs[1][i]),
			fh.write('\t'.join([list_seqs[0].split(":")[0],
			                    str(int(list_seqs[0].split(":")[1].split("-")[0])+1 + i),
			                    ".",
			                    list_seqs[1][i][1],
			                    list_seqs[2][i][ALT_BASE_NUMBER],
			                    ".",
			                    ".",
			                    ''.join([ "CTX" ,"=", ''.join([bases_in_codon[0], "(", str(list_seqs[1][i][1]), "->", list_seqs[2][i][ALT_BASE_NUMBER], ")", bases_in_codon[2] ]) ]) ,
			                    "GT", "./.\n" ]))


def parseRefGenome(RefGenFile):
	""" using RAM, and to store fasta using a dictionary """
	logger.info("Parsing Reference Genome (assuming only 5 letters possible: A,C,G,T and N) ...")
	records = list(SeqIO.parse(open(RefGenFile), "fasta"))
	dict_fasta = dict([(Rec.id, Rec.seq) for Rec in records])
	logger.info("Parsing Reference Genome (assuming only 5 letters possible: A,C,G,T and N) ... DONE")
	return dict_fasta

def parseBedFile(bedf, dict_fasta, naturalSortBedLoci=False):
	""" processing the BedFile and extract the Sequences from dico of fasta """
	lseqs = []  ## will contain all the sequences representing the targets sequences such as ATCGCCGTGTAAGCAGTGCAAGT
	lseqs_coord = []    ## will have all the coordinates of the same targets sequences listed above and in the same order as we deal with an ordered list
	logger.info("Parsing bedfile (assuming no overlapping loci) ...")
	sorted_bed_file = None
	if(naturalSortBedLoci):
		logger.info("\tSorting the Bed file Target Entries by coordinates ...")
		with open(bedf) as f:
			sorted_bed_file = natsorted(f, key=lambda y: y.lower()) ## we sort the bedfile loci in order to get the output vcf already sorted
		logger.info("\tSorting the Bed file Target Entries by coordinates ... DONE")


	for line in open(bedf) if not naturalSortBedLoci else sorted_bed_file:
		try:
			id, begin, end = line.split()   # id is the contig name, begin and end are the start and stop position respectively
			begin = int(begin) - 1 # we add one letter before the first one in order to create a trinucleotide involving the first letter in the target sequence
			end = int(end) + 1  # same comment as above but for the last base;
		except TypeError:
			raise TypeError("Expecting Integer found Something Else")
		except ValueError:
			raise ValueError("Please Check that your input BedFile is a 3-columns-tabulated file; chr start end")
		# print(':'.join([id, '-'.join([str(begin), str(end)])]))
		if id in dict_fasta:    # we check if here the contig is present in the given reference genome, if not we continue
			if (dict_fasta[id][int(begin) - 1:int(end)] is not None) and dict_fasta[id][int(begin) - 1:int(end)] != "": ## WARNING: here the -1 is not related to the coordinates, but because python starts index @ 0
				# print(dict_fasta[id][int(begin) - 1:int(end)])
				lseqs.append(str(dict_fasta[id][int(begin) - 1:int(end)]))
				lseqs_coord.append(':'.join([id, '-'.join([str(begin), str(end)])]))
	logger.debug(id + " : " + str(len(dict_fasta[id])))
	logger.info("Parsing bedfile (assuming no overlapping loci) ... DONE")
	return lseqs, lseqs_coord
#@@ End processing bedfile and getting the target coordinates and the sequences of the targets
# del sorted_bed_file
def processContigsForVcfHeader(dict_fasta):
	logger.info("Processing Contigs for vcf header ...")
	global lcontigs
	lcontigs = []  ## contains the contigs of the reference genome ;
	for id in dict_fasta:
		lcontigs.append("##contig=<ID={},length={}>".format(id, str(len(dict_fasta[id]))))
	del dict_fasta
	logger.info("\tSorting contigs ...")
	contigs = natsorted(lcontigs, key=lambda y: y.lower())
	lcontigs='\n'.join(contigs).strip()
	del contigs
	logger.info("\tSorting contigs ... DONE")
	return lcontigs
# print(lcontigs)
#
# print(len(lseqs))
# print(len(lseqs[0])) ## print the target sequences
# print(lseqs[0]) ## print the target sequences
# print(lseqs[1]) ## print the target sequences
# print(len(lseqs_coord))
# print(lseqs_coord) ## print the targets sequences normally sorted in Natural Sort

def processTargetsSequencesToCaptureAlts(lseqs, lseqs_coord):
	logger.info("processing target sequences and capturing the ALTs ...")
	pool = multiprocessing.Pool(processes=cpus)
	global pool_list_seqs
	# I know this is not clean to declare global variables, but I have been working on using Pool to make the vcfs and so far I have not found a way to make it work with multiple argument to the function wvcf()
	pool_list_seqs = None
	try:
		pool_list_seqs = pool.starmap(process_seq_for_alts, zip(lseqs, lseqs_coord))
	except Exception as err:
		logger.error(str("WARNING: Base Not Found in our set of valid Base; Base should be A, C, G, T, or N; Check Reference Genome content for non-managed bases."))
		log_err_pools.error(str("Found this base ") +str(err)+ str(" in Reference Genome; Aborting!"))
		pool.close()
		pool.terminate()
		pool.join()
		from os import kill, getpid
		kill(getpid(),9)
	finally:
		pool.close()
		pool.join()
	logger.debug(str("List of target_coordinate, codons, and Alts bases for the first target listed in user's BED file: "))
	if pool_list_seqs is not None:
		logger.debug(pool_list_seqs[0][0:3])
		logger.debug(str("Length pool_list_seqs: ") + str(len(pool_list_seqs[0])))
	logger.info("processing target sequences and capturing the ALTs ... DONE")
	return pool_list_seqs

def wvcf(N):
	with open('.m'.join([output_file_prefix, str(N+1)+".vcf"]), 'w') as of:
		of.write('\n'.join([headers, lcontigs, cmdLine, header_line]))
		write2vcf(of, pool_list_seqs, N)

def writeDataToVcfs(RefGenFile, bedf):
	logger.info("writing vcf files ...")
	global headers, cmdLine, header_line
	headers=''.join([ "##fileformat=VCFv4.1\n",
	                  "##date=\""+str(date.today().isoformat())+"\"\n",
	                  "##INFO=<ID=CTX,Number=1,Type=String,Description=\"Context (CTX) Annotation that will be use for the MutSig Coverage input file\">\n",
	                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Added this Format to make VCF perfectly to Spec and ok with Validators\">\n",
	                  "##reference="+str(RefGenFile)+"\n",
	                  "##BedFile="+str(bedf)
	                  ])
	cmdLine=''.join([ "##cmdLine=" ,"\"", str(basename(argv[0])), " ", ' '.join( argv[1:] ), "\""])
	header_line='\t'.join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SDUMMY\n"])

	## as I do not know yet if the variables that I used in the function can still be seen and well defined if I move
	## this function to the function definition sections at the top of the script, I am letting this function here for now.
	## Otherwise I need to check if we can pass several argumetn to the pool.map function
	try:
		pool = multiprocessing.Pool(processes=3)
		pool.map(wvcf, range(3))
	except Exception as err:
		logger.error("VCF writing ERROR!")
		log_err_pools.error(str("I/O error: ") + str(err))
		pool.close()
		pool.terminate()
	finally:
		pool.close()
		pool.join()
		logger.info("writing vcf files ... DONE")


if __name__ == '__main__':
	import argparse, sys, logging, multiprocessing
	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# setting up logging system
	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@
	logger = logging.getLogger('ote')
	# Define Logging parameters"
	logger.setLevel(logging.DEBUG)
	# Log to Console
	ch = logging.StreamHandler(stdout)
	# ch.setLevel(logging.INFO)
	ch.setLevel(logging.DEBUG)
	# Specifies format of log
	formatter = logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s')
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	# specifies the error level for Pools
	log_err_pools = multiprocessing.log_to_stderr()

	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# parsing and capturing user's Arguments
	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	parser = argparse.ArgumentParser(
		description='{mydoc}\n{author}, {company}, {email}'.format(mydoc=__doc__, author=__author__, company=__company__,
		                                                           email=__email__),
		formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

	pgroup = parser.add_argument_group("Input")
	# pgroup.add_argument('fasta', metavar='IN.fa|IN.fasta', type=lambda x: file_check(parser, x),
	#                     help='input1 is a reference genome file in fasta format (multi-fasta precisely)')
	pgroup.add_argument('bed', metavar='IN.bed', type=lambda x: file_check(parser, x),
	                    help='3-columns tabulatd BED file with column 1 representing the CHROMOSOME, column 2 the START and column 3 the END positions.')
	ogroup = parser.add_argument_group("Options")
	# ogroup.add_argument('-b', '--bedfile', dest='pe', action='count', default=0,
	#                     help="use paired end deduping with template. SAM/BAM alignment must contain paired end reads. Degenerate read pairs (alignments for one read of pair) will be discarded.")
	ogroup.add_argument('-f', '--fasta', dest='fasta', metavar='FASTA_FILE', type=lambda x: file_check(parser, x),
	                    default=None,
	                    help='reference genome file in fasta format (multi-fasta precisely)')
	ogroup.add_argument('-o', '--out', dest='out_prefix', default='',
	                    help='prefix of output file for sorted VCFs (default will use the input BED file as prefix and we add extension ".mX.vcf" where X will have values of 1,2 or 3 because three vcfs will be outputted by this script)')
	ogroup.add_argument('-t', '--cpus', dest='cpus', type=int, default=multiprocessing.cpu_count(), metavar='NCPUS_(int)',
	                    help="position in index read where molecular tag sequence begins. This should be a 1-based value that counts in from the 3' END of the read. (default = All cpus on the system (server|nodes|desktop))")
	ogroup.add_argument('-s', '--sort', dest='sort', action='store_true', default=False,
	                    help="By default the output vcf files will have the same contig order as the input bed file; if your bed file is not in natural sort order and you want your vcf's contigs in natural order, use this option; "
	                         "it will slow down a little bit the process; the best way is to have a bed file already sorted the way you want the vcf files sorted; it will speed up the process")
	ogroup.add_argument('-m', '--mergeLoci', dest='doLociMerge', action='store_true', default=False,
	                    help="By default we assume that there is no overlapping between the loci in the provided bed file; if there is overlapping or if not sure, use the option --mergeLoci to ensure that no overlapping exist ;"
	                         "NOTE: the output of the merging will not contain loci sorted in natural order and therefore the outpu vcfs will not be sorted in natural order as well; use the --sort to get naturally sorted vcfs")
	ogroup.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)
	# ogroup.add_argument('--errors', dest='errors', action='store_true', default=False, help=argparse.SUPPRESS)
	# ogroup.add_argument('-l', help="log file to write statistics to (optional)")
	ogroup.add_argument('-c', '--contact', action='version',
	                    version='{author}, {company}, {email}'.format(mydoc=__doc__, author=__author__, company=__company__,
	                                                                  email=__email__), help="show contact(s) and exit")
	ogroup.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
	ogroup.add_argument('-h', '--help', action='help', help='show this help message and exit')

	args = parser.parse_args()

	# Check and validate Arguments
	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	RefGenFile = args.fasta
	if args.fasta is None:
		raise ValueError(
			"Reference Genome is Missing.\nProvide A Reference Genome in Fasta format is mandatory; use -f or --fasta to provide the file as option\n")
	bedf = args.bed
	if args.out_prefix == "":
		output_file_prefix = bedf
	else:
		output_file_prefix = args.out_prefix
	cpus = args.cpus if args.cpus <= multiprocessing.cpu_count() and args.cpus > 0 else max(multiprocessing.cpu_count() - 2,
	                                                                                        1)
	naturalSortBedLoci = args.sort
	length = int(3)  # HARDCODED

	print("\n".join(["", "refGen:\t" + RefGenFile, "bedfile:\t" + bedf, "prefix:\t" + output_file_prefix, "cpus:\t" + str(cpus),
	                 "sort:\t" + str(naturalSortBedLoci), ""]))

	try:
		if (args.doLociMerge):
			bedf = mergeLoci(bedf)
		dict_fasta = parseRefGenome(RefGenFile)
		lseqs, lseqs_coord = parseBedFile(bedf, dict_fasta, naturalSortBedLoci)
		lcontigs = processContigsForVcfHeader(dict_fasta)
		pool_list_seqs = processTargetsSequencesToCaptureAlts(lseqs, lseqs_coord)
		writeDataToVcfs(RefGenFile, bedf)
		logger.info("Job completed successfully")
	except Exception as e:
		logger.info(e.message)
		exit(1)
	exit()
