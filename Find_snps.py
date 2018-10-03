#!/usr/bin/python3
##second step in snp discovery from GT-seq reads pipeline
##Looks through haplotypes and finds potential biallelic markers 
#all three input parameters are optional
#default -min_prop is 0.8
#default -min_maf is 0.05
#default -min_depth is 10
#In terminal type python Find_snps.py -min_prop min_proportion_of_samples_with_haplotypes -min_maf min_minor_allele_frequency -min_depth min_depth_for_calling_haplotypes


from math import ceil
import sys

def Main():

	samples = {}
	#read in sample names, marker names, and haplotypes
	#all_haplotypes.txt is tab-delimited indiv_name, marker_name, hap_1, depth_1, hap_2, depth_2
	for line in open('all_haplotypes.txt', 'r'):
		line = line.rstrip()
		sep = line.split('\t')
		if sep[0] in samples:
			samples[sep[0]][sep[1]] = sep[2:6]
		else:
			samples[sep[0]] = {sep[1] : sep[2:6]}

			
	markers = []		#get list of markers from one individual (all individuals will have entries for all markers, whether they were typed for them or not)
	for i in samples[list(samples.keys())[0]]:
		markers.append(i)
	
	min_num_haplo = ceil(0.8*len(samples))*2		#use default of 0.8 if no value is given
	min_maf = 0.05		#use default of 0.05 if no value is given
	min_depth = 10		#use default of 10 if no value is given
	
	for flag in range(1, len(sys.argv), 1):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == '-min_depth':
			min_depth = int(sys.argv[flag + 1])
		if sys.argv[flag] == '-min_maf':
			min_maf = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-min_prop':
			min_num_haplo = ceil(float(sys.argv[flag + 1])*len(samples))*2

	print('Finding markers with reads in at least ', str(int(min_num_haplo/2)), ' samples\n\t a minimum minor allele frequency of ', min_maf, '\n\t and using a minimum depth of ', min_depth)
	
	poten_out = open('potential_snps.txt', 'w')
	poten_out.write('locus\tcol\tmaf\n')
	for marker in markers:	#for all markers
		#get all haplotypes for the marker
		all_haplotypes = []
		for indiv in samples:	#and for all samples
			if samples[indiv][marker][0][0] != 'N' and int(samples[indiv][marker][1]) >= min_depth: #check that there is an actual read and that depth is above min_depth
				all_haplotypes.append(samples[indiv][marker][0])		
			if samples[indiv][marker][2][0] != 'N' and int(samples[indiv][marker][3]) >= min_depth: 
				all_haplotypes.append(samples[indiv][marker][2])
				
		#discard markers that didn't have haplotypes in at least x% of individuals
		if len(all_haplotypes) <= min_num_haplo:
			continue
		#find snps 
		max_length = 0
		for i in all_haplotypes:
			if len(i) > max_length:
				max_length = len(i)		#find maximum length of haplotypes
		for i in range(0, max_length, 1):
			unique = {}					#build dictionary of alleles and their frequency at a given position in the haplotype
			for haplo in all_haplotypes:
				if len(haplo) < (i + 1):		#make sure current haplotype is long enough to have a base
					continue
				if haplo[i] == 'N':		#prevent no call from being counted as an allele
					continue
				if haplo[i] in unique:
					unique[haplo[i]] += 1
				else:
					unique[haplo[i]] = 1
			if len(unique) != 2:		#check that locus is biallelic (not monomorphic and not more than 2 alleles)
				continue
			#calculate maf
			first = True
			for j in unique:
				if first:
					minor = j
					first=False
				elif unique[j] < unique[minor]:
					major = minor
					minor = j
				else:
					major = j
			maf = unique[minor]/(unique[major] + unique[minor])
			
			#output snp and maf if it passes maf cutoff
			if maf >= min_maf:
				poten_out.write(marker + '\t' + str(i) + '\t' + str(maf) + '\n')
	
	poten_out.close()

Main()


