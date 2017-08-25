# CovGen
Creates a capture specific exome_full192.coverage.txt file required by MutSig

## Script Purpose 
Create Three VCF files with alternates mutations to each base of the target sequences in bed File ; 
These VCF files  will further be annotated by snpEff (or any vcf annotator) to get the Effect and 
use that effect to create the coverage file for MutSig  

#### PreRequisites to speed up the process:
* Loci in BED File should be sorted the same way you would like to have the position sorted in the 3 output vcf files; if not use the `--sort` and the script will sort the positions in natural order
* The input BED file __cannot__ have overlapping loci; If this is the case, duplcaited lines in output vcfs will be createdand will cause issues downstream analysis; if you know that the vcf has overlapping positions, use `--merge` options to the script.
***
## System requirements
* python 3.4 or above
* modules:
  + [pybedtools 0.7.7 or up](https://pypi.python.org/pypi/pybedtools/0.7.9)
  + multiprocessing
  + [BioPython](http://biopython.org/wiki/Download)
  + [natsort](https://pypi.python.org/pypi/natsort)
  + sys
  + os
  + logging
  + functools



***
## Usage
To get the arguments and options to the script, run:  

` python getALTsForTargetsSeqsForMutSig.py --help `  or  `python getALTsForTargetsSeqsForMutSig.py -h`


minimum requirements are to be provided to the script:
* reference genome (GRCh38 or older)
* bed file containing target positions related to the provided reference genome; __3-columns tabulated bed file__
***

www.tgen.org

