#!/usr/bin/env python
##first step in snp discovery from GT-seq reads pipeline
##requires fastq files form all individuals to be in one folder and to be named: individualname.fastq
##requires giving the file path to the probeseq in the terminal
##optional multiprocessing - give the number of threads you want to use (default is teh number of threads the machine has)
## differentiates het's and hom's by ratio of read counts, (top depth) / (second top depth) with hets having the ratio less than or equal to het_ratio (default is 10)
##		specify het_ratio using the -r flag on the command line
#sorts reads in each individual fastq by fwd primer (from probeseq), 
#determines top two sequences (haplotypes) with highest depth, outputs both if ratio of depths is < het_ratio (heterozygote), else outputs only top sequence (homozygote),
#writes a tab-sep-value file "all_haplotypes.txt" containing individualname	locusname	haplotype_1	depth_1	haplotype_2	depth_2
#in terminal type: python Sort_reads.py -ps file/path/to_probeseq.csv -t number_of_threads -r het_ratio
#
### this copy does not use the biopython library
#

import re, glob, sys
from multiprocessing import Process, Queue, cpu_count

def Sort(queue, samples, name_fwd, het_ratio):
	haplo_process = []
	for indiv in samples:
		haplo_indiv = []
		for marker in name_fwd:		#iterate through markers - could assess all markers at the same time, but memory demand might be too high
			length = len(name_fwd[marker])	#get length of fwd primer
			unique = {}			#dictionary to count depth of each sequence
			with open(indiv + '.fastq', 'r') as fastq_file:
				next(fastq_file) #skip header
				string_read = fastq_file.readline().rstrip()#read sequence
				while string_read:
					if name_fwd[marker] == string_read[0:length]:
						if string_read in unique:
							unique[string_read] += 1	#add one
						else:
							unique[string_read] = 1	#add new transcript to dict
					next(fastq_file) # skip spacer
					next(fastq_file) #skip  quality scores
					string_read = fastq_file.readline()	#read header to determine if there is another one
					string_read = fastq_file.readline().rstrip()	#read next seq
								
			#output top two sequences if ratio is less than 10, top one if ratio is higher than 10
			first = ['N', 0.01]
			second = ['N', 0.01]		#0.01 to prevent division by zero if no reads
			for sequence in unique:
				if unique[sequence] > first[1]:
					second = first
					first = [sequence, unique[sequence]]
				elif unique[sequence] > second[1]:
					second = [sequence, unique[sequence]]
			if first[1]/second[1] <= het_ratio:							#heterozygous
				haplo_indiv.append([marker, first[0], first[1], second[0], second[1]])
			else:
				haplo_indiv.append([marker, first[0], first[1], first[0], first[1]])	#if homozygous, write same sequence twice, in order to make reading 
																						#in genotypes later easier (no need to test length to determine if indiv is het or homo)
		for i in haplo_indiv:		#add haplotypes to running list for the process
			haplo_process.append([indiv] + i)
	queue.put(haplo_process)

def Main():
	#read in list of samples
	samples = glob.glob('*.fastq')
	samples = [re.sub('\.fastq', '', x) for x in samples]
	num_samples = len(samples)
	
	
	threads = cpu_count()	#default values if no value is given
	ps_file = 'no_input'
	ratio = 10
		
	for flag in range(1, len(sys.argv), 1):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == '-t':
			threads = int(sys.argv[flag + 1])
		if sys.argv[flag] == '-ps':
			ps_file = sys.argv[flag + 1]
		if sys.argv[flag] == '-r':
			ratio = float(sys.argv[flag + 1])

	
	#read in probeseq
	if ps_file == 'no_input':
		print('No probeseq file specified. Exiting.')
		return
		
	markers = {}
	try:
		for line in open(ps_file, 'r'):
			separated = line.split(",")
			markers[separated[0]] = separated[5]
	except:
		print('Could not open probeseq specified. Exiting')
		return
		
	
	#divide samples up amongst number of processes to use
	
	print('Using ', str(threads), 'threads')
	per_process = num_samples//threads
	leftover = num_samples % threads
	thread_dict = {}
	count = 0
	for i in range(0, threads, 1):
		if i < leftover:
			thread_dict[i] = samples[count : count + per_process + 1]
			count += per_process + 1
		else:
			thread_dict[i] = samples[count : count + per_process]
			count += per_process
	
	#start processes of sorting
	processes = {}
	queues = {}
	for i in thread_dict:
		queues[i] = Queue()
		processes[i] = Process(target=Sort, args=(queues[i], thread_dict[i], markers, ratio))
		processes[i].start()
	
	#get returned haplotypes and depths
	with open('all_haplotypes.txt', 'w') as out_file:
		for i in queues:
			results = queues[i].get()
			for j in results:						#write returned haplotypes to a single file for use by later scripts
				out_file.write(j[0] + '\t' + j[1] + '\t' + j[2] + '\t' + str(j[3]) + '\t' + j[4] + '\t' + str(j[5]) + '\n') 
	
	#join processes
	for i in processes:
		processes[i].join()

Main()