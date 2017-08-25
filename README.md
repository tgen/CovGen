# CovGen
Creates a capture specific exome_full192.coverage.txt file required by MutSig

## Summary
MutSig provides a "territory" table (exome_full192.coverage.txt) for times when detailed coverage information is not available for each sample in your cohort. This coverage file may not properly represent the target space utilized by your capture kit and can adversely affect the results of your mutsig analysis.
CovGen bridges the gap between detailed sample level coverage information and the exome_full192.coverage.txt table that MutSig provides with a target specific full coverage table.

There are a few fundamental differences or caveats between the coverage provided by MutSig and the one produced by CovGen. ENSG Ensembl IDs are used in place of HUGO symbols. This requires that the covariates file utilized by MutSig must also be converted to Ensembl IDs. After MutSig analysis the Ensambl ID's can easily be mapped back to HUGO ID's for readability. CovGen only considers protein coding genes as defined by CDC feature type in the user provided Ensembl GTF. Alternate alleles that are upstream or downstream for a given gene are excluded. 

In addition to the coverage file, CovGen also outputs a bed file representing the final target space used to create the coverage file and an ENSG list. These two files should be used to filter your mutation file (MAF). This step helps to prevent MutSig from passing the following warning and zeroing out all noncoding mutations and coverage for the rest of the calculation. 

`WARNING:  coding and noncoding rates are too different`

If you have annotated your variants using snpEff with the ANN annotation standard then the snpEff_ANN_mutation_type_dictionary_file.txt provided in this package can be used in place of the mutation_type_dictionary_file.txt provided by MutSig. It is a good idea to review the mapping of the ANN Variant_Classification to MutSig effects as a few of the mappings could be open to interpretation.

## System requirements
* 36G of ram
* 8 processing cores
* [samtools 1.5 or above](http://www.htslib.org/download/)
* [bedtools 2.25.0 or above](http://bedtools.readthedocs.io/en/latest/content/installation.html)
* [snpEff 4.2 or above](http://snpeff.sourceforge.net/)
* [gawk 3.1.7 or above](https://www.gnu.org/software/gawk/)
* python 3.4 or above with the following modules:
    + [pybedtools 0.7.7 or up](https://pypi.python.org/pypi/pybedtools/0.7.9)
    + multiprocessing
    + [BioPython](http://biopython.org/wiki/Download)
    + [natsort](https://pypi.python.org/pypi/natsort)
    + sys
    + os
    + logging
    + functools

***
## Options

| Option  | Argument  | Required  | Description |
| ------- |:--------- |:---------:|:-------------- |
| -o      | string  |Yes|  Prefix for output files|
| -f      |file     |Yes| Reference genome fasta file |
| -g      |file     |Yes| Ensembl annotations in GTF format |       
| -t      | file    |Yes|  Non-Padded zero based exome capture targets bed file. CovGen pads each target by 100bp on each end for you. |      
| -s | path      |Yes| Full path to snpEff directory. Other annotators are not currently supported.     |
| -v | string   |Yes| SnpEff genome_version to use. Other annotators are not currently supported.     |     
| -b | file      |No| File of full paths to 6 BAM files. Some probes are not as effective as others. This option uses 6 bam files to filter out bases from your targets when 2 or more samples have <=10x coverage. For chrY, 5 samples with <=10x get filtered out. Please include at least 2 female samples to prevent all of Y being excluded.  |      
| -p | int       |No| Number of processors available. The VCF generation step requires the number of processors available to be set. All other steps will automatically use all cores available to the processor. Default = 3   |
| -e | file      |No| List of ENSGs to filter out.                 |

## Usage
To get the arguments and options to the script, run:  

`CovGen `  or  `CovGen --help`

Basic usage with required options

```
CovGen -o Agilent_SureSelect_V5_plusUTR \
  -f hs37d5.fa \
  -g Homo_sapiens.GRCh37.74.gtf \
  -t zeroBased_targets.bed \
  -s /path/to/snpEff/directory \
  -v GRCh37.74 
```
Usage with additional filtering options and number of processes for vcf generation tool

```
CovGen -o Agilent_SureSelect_V5_plusUTR \
  -f hs37d5.fa \
  -g Homo_sapiens.GRCh37.74.gtf \
  -t zeroBased_targets.bed \
  -s /path/to/snpEff/directory \
  -v GRCh37.74 \
  -p 25 \
  -b list_of_six_bam_files.txt \
  -e ENSG_list_to_filter_out.txt 
```


www.tgen.org

