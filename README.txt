#############################
Discovering new SNPs from existing GT-seq panels

#############################

This series of scripts was written to be used in the scenario where you already have a GT-seq primer mix for one species, but you are interested
in either genotyping a population that isn't polymorphic (or has low maf's) for the SNPs you are typing usign that primer mix, or you are
genotyping a closely related species, and so you get reasonable amplification with the primer mix but not the SNPs that are being
typed are monomorphic. In these scenarios, if you could find other SNPs in the loci that are being amplified by your primer mix, you could avoid
going through the entire process of panel development. This series of three scripts is aimed at finding these SNPs.

First, you need to run GT-seq with the primer mix and a group of the samples from the population you want to find SNPs in. Separate out the reads into
separate FASTQ files for each individual, named samplename.fastq, and placed in one directory. Save all scripts to this directory as well.

The first script (Sort_reads.py) will sort through the reads for each individual and group all the reads by the forward primer (from the probeseq file you specify)
that they contain. It then finds the two unique sequences with highest read count in each individual/forward primer combination. If the ratio of read counts between these (higher read count / lower read count) is greater 
than the specified het_ratio (default is 10), it considers the individual homozygous for the more common haplotype, if less or equal, it considers the individual heterozygous. All haplotypes for
all loci for all individuals are output in "all_haplotypes.txt", which is used by the second and third scripts.

The second script (Find_snps.py) compares all of the haplotypes within each locus, and determines where there are differences (considered as potential SNPs) within the 
group of samples. It outputs any differences and the minor allele frequency at that locus in a file called "potential_snps.txt". This file contains the locus name (locus), 
the location (col) of the potential SNP in the locus (0-based numbering, first base is 0) and the minor allele frequency (maf) of that potential SNP. This file is then used 
by the third script.

The third script (Build_probeseq.py) builds a probeseq file for you to use to test the efficacy of the potential SNPs. It combs through the potential_snps.txt file and identifies one 
SNP for each locus to put into the probeseq. It chooses the SNP with the highest MAF, ties being broken by choosing the SNP closer to the forward primer. The 
only exception to this is that if you have a file named "skip_snps.txt" in the directory that contains SNPs in the format locusname_col it will skip these SNPs. 
It will then create the in silico probes for this locus, using variation at other sites in the probe as variable bases in the in silico probes. It then outputs a probeseq 
file for you to use in testing these SNPs. 

The scripts are best at finding substitutions. Indels will be found by the second script, but they will show up as a long run of adjacent SNPs, and not as one indel. Additionally, 
they will most likely be filtered out by the third script, unless you make the -h parameter close to 1 (but you may also encounter other problems from this).

After you generate a probeseq, it is important to test the markers by analyzing your GT-seq reads with the GT-seq pipeline and the new probeseq file. Some things to watch 
out for are: 1. SNPs that turn out to not be polymorphic. The scripts are written to be more liberal when deciding if something is a potential SNP or not, with the justification that a 
false positive will be filtered out quickly, but a false negative will be a loss of a marker; 2. SNPs for which all, or almost all, individuals appear to be heterozygous.
This can occur when two different locations in the genome are being amplified by the same forward primer, and the "SNP" that was discovered is simply a difference between these
two locations of the genome. In both of these cases, the "SNPs" can easily be deleted from the probeseq, either manually, or using the optional skip_snps.txt file. Once you have
decided on a group of SNPs to use, you may also need to adjust allele correction values for the new SNPs to optimize your genotyping success. 

Specific explanations of all the options for the three scripts are below.

#############################

python Sort_reads.py -ps file/path/to_probeseq.csv -t number_of_threads -r het_ratio
	-ps give the file path to the probeseq file for the GT-seq panel you are using to find new SNPs. The script will use the locus name and the forward
			primer in the probeseq file.
	-t give the number of threads you would like the program to use, if you want to use multiple threads (default: 1)
	-r give the ratio of read counts used to determine if a sample is homozygous or heterozygous 


python Find_snps.py -min_prop min_proportion_of_samples_with_haplotypes -min_maf min_minor_allele_frequency -min_depth min_depth_for_calling_haplotypes
	-min_prop give the minimum proportion of samples with haplotypes necessary to consider potential SNPs within a locus (default: 0.8)
	-min_maf give the minimum minor allele frequency to consider a potential SNP (default: 0.05)
	-min_depth give the minimum number of reads to call a haplotype (default: 10)


python Build_probeseq.py -v minimum_prop_of_haplotypes_for_variable_bases_in_probe -ps file/path/to/probeseq -min_depth min_depth_to_call_haplotype -h max_prop_of_seq_with_snps -suf suffix_for_output_probeseq
	-v give the minimum frequency of an alternative base in the probe sequence (other than the SNP of interest) to include that base as a variable base in the probe (default: 0)
	-min_depth give the minimum number of reads to call a haplotype (default: 10)
	-h give the proportion to use in calculating the maximum number of potential SNPs in a given locus for SNPs in the locus to be considered. The maximum number of SNPs is 
		calculated as (proportion_from_-h * max_read_length_for_locus_in_your_data). This parameter is helpful for removing loci from forward primers that are amplifying two 
		locations in the genome (default: 0.15)
	-suf give the suffix to append to the filename of the output probeseq. Will be named "probeseq_suffix.csv" (default: no suffix)

	
skip_snps.txt file format:

locus1_col1
locus1_col2
locus5_col1

for example:

Ots_23577-43_78
Ots_26691-36_65
Ots_30619-61_59
Ots_3209-10_78
Ots_37492-53_74
