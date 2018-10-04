#!/usr/bin/python3
##third step in snp discovery from GT-seq reads pipeline
##Looks through haplotypes and potential_snps.txt to build a probeseq file
##will look for a file named skip_snps.txt to skip over previously attempted snps
##skip_snps.txt should be list of snps to skip with each snp on a new line and written as "locus_col" with locus and col matching the potential_snps.txt file
# -v default is 0
# -min_depth default is 10
# -h default is 0.15
# -suf default is none
#in terminal type python Build_probeseq.py -v minimum_prop_of_haplotypes_for_variable_bases_in_probe -ps file/path/to/probeseq -min_depth min_depth_to_call_haplotype -h max_prop_of_seq_with_snps -suf suffix_for_output_probeseq

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

	poten = {}
	try:
		first = True
		for line in open('potential_snps.txt', 'r'):
			if first:
				first = False
				continue
			sep = line.rstrip('\n').split('\t')
			if sep[0] in poten:
				poten[sep[0]].append([int(sep[1]), float(sep[2])])
			else:
				poten[sep[0]] = [[int(sep[1]), float(sep[2])]]			#create dict "poten" with key as marker name, value as array  [[snp1_position_0_based, maf], [snp2_position_0_based, maf]]
	except:
		print('cannot find/read potential_snps.txt. Exiting')
		return
	
	vari = 0	#default values if no value is given
	min_depth = 10
	homol = 0.15
	suf = ''
	
	for flag in range(1, len(sys.argv), 1):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == '-min_depth':
			min_depth = int(sys.argv[flag + 1])
		if sys.argv[flag] == '-v':
			vari = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-ps':
			ps_file = sys.argv[flag + 1]
		if sys.argv[flag] == '-h':
			homol = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-suf':
			suf = str(sys.argv[flag + 1])
	
	#read in probeseq
	probeseq = {}
	try:
		for line in open(ps_file, 'r'):
			separated = line.split(',')
			probeseq[separated[0]] = separated[5]
	except:
		print('Could not open probeseq specified. Exiting')
		return	
	
	#read in skip file, if present
	skip = []
	try:
		for line in open('skip_snps.txt', 'r'):
			skip.append(line.rstrip())
		print('Read in ', str(len(skip)), ' snps to skip.')
	except:
		print('Did not find skip_snps.txt. Will choose all snps based on highest maf.')
		
	print('Incorporating variation present at a frequency equal to or greater than', vari, 'into the probes')
	print('Using a minimum read depth of', min_depth, 'to call haplotypes')
	
	#choose snps based on highest maf
	chosen = {}
	for i in poten:
		temp = [0,0]
		for j in poten[i]:
			if temp[1] < j[1] and str(i) + '_' + str(j[0]) not in skip:	#choose highest maf that isn't skipped
				temp = j
		if temp == [0,0]:	#if all snps were skipped, continue to the next marker
			print('All potential snps in marker ', i, 'were skipped')
			continue
		chosen[i] = temp

	#build probeseq - Marker_name, allele1, allele2, probe1, probe2, fwd, corr1, corr2
	probeseq_out = open('probeseq_' + suf + '.csv', 'w')
	for locus in chosen:
		#get all haplotypes for the marker
		all_haplotypes = []
		for indiv in samples:	#for all samples
			if samples[indiv][locus][0][0] != 'N' and int(samples[indiv][locus][1]) >= min_depth: #check that there is an actual read and that depth is above min_depth
				all_haplotypes.append(samples[indiv][locus][0])		
			if samples[indiv][locus][2][0] != 'N' and int(samples[indiv][locus][3]) >= min_depth: 
				all_haplotypes.append(samples[indiv][locus][2])
		max_length = 0
		for i in all_haplotypes:
			if len(i) > max_length:
				max_length = len(i)		#find maximum length of haplotypes
		if len(poten[locus]) >= homol * max_length:			#check if potential amplification of more than one locus - note: this will also remove indels
			print('Skipping locus', locus, 'due to potentailly amplifying multiple regions')
			continue
		vari2 = vari * len(all_haplotypes)	#expand 
		#build probes
		snp_pos = chosen[locus][0]
		if snp_pos + 8 <= max_length:
			start = snp_pos - 6
			end = snp_pos + 8
		else:
			end = max_length
			start = max_length - 14
		probe1, probe2 = '', ''
		#get all alleles and their frequency throughout probes
		for i in range(start, end, 1):
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
			if i == snp_pos:			#give different alleles to each probe
				first = True
				for allele in unique:
					if first:
						probe1 += allele
						A1 = allele			#save allele for output to probeseq
						first = False
					else:
						probe2 += allele
						A2 = allele			#save allele for output to probeseq
			else:
				temp =[]
				for allele in unique:
					if unique[allele] >= vari2:	#add variation greater than vari, excluding the snp of interest
						temp.append(allele)
				if len(temp) == 1:
					probe1 += temp[0]
					probe2 += temp[0]
				else:
					probe1 += '['
					probe2 += '['
					for alt in temp:
						probe1 += alt
						probe2 += alt
					probe1 += ']'
					probe2 += ']'
		#turn probes into strings
		
		#output line of probeseq
		probeseq_out.write(locus + '_' + str(snp_pos) + ',' + A1 + ',' + A2 + ',' + probe1 + ',' + probe2 + ',' + probeseq[locus] + ',' + '0' + ',' + '0' + '\n')
		
	probeseq_out.close()
	
Main()
