#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "19/05/2015"
__version__ = "$Revision: 2.5"

from optparse import OptionParser
from operator import itemgetter
from collections import Counter, deque
from itertools import chain
from multiprocessing import Pool, Queue, Process
import sys, time, os
import random
import copy
import json, csv
import cPickle as pickle
import shelve
import resource
import ntpath
import tarfile

""" global variables """
# LSH
leftPartRatio = 0.5
rightPartRatio = 0.5
printMinhashProcess = 5000000
setVersion = False

# Sequence Alignment
M1 = 1
M2 = 2
secondSample = 0
overlap = 6 # Overlap region in both directions i.e. 20 overlap in total if 10
maxAlignments = 1 # per read
requiredOverlaps = 3
maxCandMates = 5000
MUTFIND = 1
p_id = -1


# Test variables
c1 = 0
c2 = 0
c3 = 0
c4 = 0
c5 = 0
c6 = 0
c7 = 0
c8 = 0
numreadL = 0
numreadR = 0

# ************************************************************************** #
#																			 #
#									Helpers									 #
#																			 #
# ************************************************************************** #
def optionParse():
	"""
	Parse arguments from command line.
	"""
	desc = """Compare sets of sequencing data to find mutations."""

	parser = OptionParser(usage="usage: %prog --fasta_file filename",
						  description=desc,
						  version="%prog version 2.0")

	parser.add_option("-f", "--fasta_file",
					  metavar="<FILENAME>",
					  default="../Data/Fasta/reads.fa",
					  action="store",
					  dest="fasta_file",
					  help="set <FILENAME> as fasta file.")

	parser.add_option("-n", "--normal",
					  metavar="<FILENAME>",
					  action="store",
					  dest="normal_file",
					  help="set <FILENAME> as normal reads sample.")

	parser.add_option("-d", "--diseased",
					  metavar="<FILENAME>",
					  action="store",
					  dest="diseased_file",
					  help="set <FILENAME> as diseased reads sample.")

	parser.add_option("-i", "--candidate_pairs",
					  metavar="<FILENAME>",
					  type=str,
					  action="store",
					  dest="input",
					  help="set <FILENAME> as input file containing \
							candidate pairs to import.")

	parser.add_option("-e", "--export_candidate_pairs",
					  metavar="<FILENAME>",
					  type=str,
					  action="store",
					  dest="output",
					  help="set <FILENAME> as output file to export \
							candidate pairs to. Dump as either txt, json, \
							pickle or csv file. Just put the right extension\
							and everything else is automated.")

	parser.add_option("-l", "--log_file",
					  metavar="<FILENAME>",
					  type=str,
					  default="log.txt",
					  action="store",
					  dest="log",
					  help="set <FILENAME> as log file.")

	parser.add_option("-k", "--k_size",
					  metavar="<VALUE>",
					  type=int,
					  default=12,
					  action="store",
					  dest="k",
					  help="set <VALUE> as size for k-shingles.")

	parser.add_option("-t", "--threshold",
					  metavar="<VALUE>",
					  type=float,
					  default=float(2)/3,
					  action="store",
					  dest="threshold",
					  help="set <VALUE> as threshold similarity.")

	parser.add_option("-b", "--bands",
					  metavar="<Value>",
					  type=int,
					  default=25,
					  action="store",
					  dest="bands",
					  help="set <VALUE> as the number of bands for LSH.")

	parser.add_option("-r", "--rows",
					  metavar="<Value>",
					  type=int,
					  default=40,
					  action="store",
					  dest="rows",
					  help="set <VALUE> as the number of rows for LSH.")

	parser.add_option("-x", "--similarity_measure",
					  metavar="<VALUE>",
					  type=str,
					  default="naive",
					  action="store",
					  dest="sim",
					  help="set <VALUE> as similairy measure to use (obsolete).")

	parser.add_option("-m", "--minhash_alg",
					  metavar="<VALUE>",
					  type=int,
					  default="6",
					  action="store",
					  dest="m",
					  help="<VALUE> defines the minhash algorithm to use.")

	parser.add_option("-s", "--seed",
					  metavar="<VALUE>",
					  type=int,
					  default=42,
					  action="store",
					  dest="s",
					  help="set <VALUE> as seed for hash functions.")

	parser.add_option("-T", "--test",
					  metavar="<VALUE>",
					  type=int,
					  action="store",
					  dest="T",
					  help="perform test <VALUE>.")

	parser.add_option("-S", "--supporting_reads",
					  metavar="<VALUE>",
					  type=int,
					  default=3,
					  action="store",
					  dest="S",
					  help="Reads required to support mutation.")

	parser.add_option("-o", "--mismatch_overlap",
					  metavar="<VALUE>",
					  type=int,
					  default=1,
					  action="store",
					  dest="o",
					  help="Number of allowed mismatches in overlap.")

	parser.add_option("-g", "--mismatch_group",
					  metavar="<VALUE>",
					  type=int,
					  default=2,
					  action="store",
					  dest="g",
					  help="Number of allowed mismatches in group.")



	(options, args) = parser.parse_args()

	global requiredOverlaps
	requiredOverlaps = options.S
	global M1
	M1 = options.o
	global M2
	M2 = options.g

	return options.fasta_file, options.normal_file, options.diseased_file,\
		   options.k, options.threshold, options.bands, options.rows,\
		   options.m, options.s, options.log, options.input, options.output,\
		   options.T


def memory_usage_resource():
	"""
	Computes the ressource usage (in MB) at a given time during runtime.
	Computes total amount used so far, so not the amount currently in use.
	"""
	rusage_denom = 1024.
	if sys.platform == 'darwin':
		# ... it seems that in OSX the output is different units ...
		rusage_denom = rusage_denom * rusage_denom
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
	return mem


def logprint(log_file, flush, *output):
	"""
	Prints to standard out and to log file
	"""
	log_output = ""
	log_output = " ".join(map(str, output))
	print log_output
	if log_file:
		log_file.write(log_output + "\n")
		if flush:
			log_file.flush()


def exportCandidatePairs(candidatePairs, output_file, log, numReads=None):
	"""
	Export candidate pairs to a file.
	The type of file is determined on the provided filename for output_file.
	Supported filetypes: txt, json, pickle (python) and csv.
	"""
	tim = time.clock()
	# Output file extension
	ext = output_file.rsplit(".", 1)
	# default to txt if no extension provided
	if len(ext) == 0:
		ext = "txt"
		output_file += ".txt"
	else:
		ext = ext[1]

	# save set information - However, more space consuming
	# and not needed. Hence, this should never be used.
	if ext == "set_pickle":
		with open(output_file, "w") as f:
			pickle.dump(candidatePairs, f)

	elif ext == "json":
		with open(output_file, "w") as f:
			if isinstance(candidatePairs[0], set):
				for id1 in candidatePairs:
					candidatePairs[id1] = list(candidatePairs[id1])
			json.dump(candidatePairs, f)

	elif ext == "pickle":
		with open(output_file, "w") as f:
			if isinstance(candidatePairs[0], set):
				for id1 in candidatePairs:
					candidatePairs[id1] = list(candidatePairs[id1])
			pickle.dump(candidatePairs, f)

	elif ext == "txt":
		with open(output_file, "w") as f:
			for id1 in candidatePairs:
				f.write(str(id1)+"\t")
				#sortedElements = sorted(list(candidatePairs[id1]))
				sortedElements = list(candidatePairs[id1])
				for id2 in sortedElements[:-1]:
					f.write(str(id2)+",")
				if len(sortedElements) > 0:
					f.write(str(sortedElements[-1])+"\n")

	# Test-only write the first numReads reads to output file.
	elif ext == "temp":
		with open(output_file, "w") as f:
			for id1 in xrange(numReads):
				f.write(str(id1)+"\t")
				#sortedElements = sorted(list(candidatePairs[id1]))
				#print sortedElements
				if id1 in candidatePairs:
					sortedElements = list(candidatePairs[id1])
					for id2 in sortedElements[:-1]:
						f.write(str(id2)+",")
					if len(sortedElements) > 0:
						f.write(str(sortedElements[-1]))
				f.write("\n")

	elif ext == "txt2":
		with open(output_file, "w") as f:
			for id1 in candidatePairs:
				for id2 in candidatePairs[id1]:
					f.write(str(id1)+"\t"+str(id2)+"\n")

	elif ext == "csv":
		w = csv.writer(open(output_file+".csv", "w"))
		for key, val in candidatePairs.items():
			w.writerow([key, val])

	# Else export to whatever filename that is provided in the format
	# used for txt files.
	else:
		output_file += ".txt"
		with open(output_file, "w") as f:
			for id1 in candidatePairs:
				f.write(str(id1)+"\t")
				sortedElements = list(candidatePairs[id1])
				for id2 in sortedElements[:-1]:
					f.write(str(id2)+",")
				f.write(str(sortedElements[-1])+"\n")

	logprint(log, False, "Exported candidate pairs to", output_file,
			 "in", time.clock()-tim, "seconds")


def importCandidatePairs(input_file, log):
	"""
	Import candidate pairs from a file.
	Supported filetypes: txt, json, pickle (python) and csv.
	"""
	tim = time.clock()
	logprint(log, True, "Importing candidate pairs...")
	candidatePairs = dict()
	# Input file extension
	ext = input_file.rsplit(".", 1)
	if len(ext) == 0:
		ext = "txt"
		input_file += ".txt"
	else:
		ext = ext[1]

	if ext == "set_pickle":
		with open(input_file, "r") as f:
			candidatePairs = pickle.load(f)
		logprint(log, False, "pickle-set:", time.clock()-tim)

	elif ext == "json":
		with open(input_file, "r") as f:
			candidatePairs = json.load(f)
		logprint(log, False, "json:", time.clock()-tim)

	elif ext == "pickle":
		with open(input_file, "r") as f:
			candidatePairs = pickle.load(f)
		logprint(log, False, "pickle-list:", time.clock()-tim)

	elif ext == "txt":
		with open(input_file, "r") as f:
			imported = 0
			for line in f:
				elements = line.split()
				key = int(elements[0])
				pairs = map(int, elements[1].split(','))
				candidatePairs[key] = pairs
				imported += len(pairs)
				# print elements
				# print key
				# print pairs
				# print candidatePairs[key]
				# sys.exit()
		logprint(log, False, "Imported", imported/2, "candidate pairs from",
				 input_file, "in", time.clock()-tim, "seconds.")

	elif ext == "txt2":
		with open(input_file, "r") as f:
			for line in f:
				elements = map(int, line.split())
				if elements[0] in candidatePairs:
					candidatePairs[elements[0]].append(elements[1])
				else:
					candidatePairs[elements[0]] = [elements[1]]
		logprint(log, False, "txt2 file:", time.clock()-tim)

	elif ext == "csv":
		for key, val in csv.reader(open(input_file)):
			if key in candidatePairs:
				candidatePairs[key].append(val)
			else:
				candidatePairs[key] = [val]
		logprint(log, False, "csv:", time.clock()-tim)

	else:
		logprint(log, True, "File format is not supported for input file."
				 "Please specify file format (extension) as either txt,",
				 "json, pickle or csv.")
		sys.exit()

	# Print imported candidate pairs
	# for id1 in candidatePairs:
	#	  print id1
	#	  print candidatePairs[id1]
	#	  sys.exit()

	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	return candidatePairs


# ************************************************************************** #
#																			 #
#							   Pre-computations								 #
#																			 #
# ************************************************************************** #
def computeHashFunctions(n, shingles, log):
	"""
	Computes n lists of shuffled elements from 1..#shingles.
	These lists represents the hash functions needed for LSH.
	"""
	# Create n different permutations (hash functions) of the shingles
	tim = time.clock()
	hashfuncs = []
	for i in xrange(n):
		h = range(len(shingles))
		random.shuffle(h)
		hashfuncs.append(h)
		# print h,"\n"
	logprint(log, False, "Computed hashfunctions in",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	return hashfuncs


def computeShinglesTable(fasta_file, shinglesPos, k, log):
	"""
	Computes a table for fast look-up of k-shingles and their corresponding
	position in the reads - when it was first encountered in the fasta file.
	"""
	if not fasta_file:
		return shinglesPos
	logprint(log, True, "Computing table of shingles positions...")
	tim = time.clock()
	with open(fasta_file, "rU") as fasta_file:
		#shinglesPos = dict()
		read = ""
		pos = 0
		for line in fasta_file:
			# If line starts with ">", which indicates end of a sequence
			# Then, append it to list of reads
			if line.startswith(">"):
				if read != "":
					# Splits the string into two parts
					leftpart = read[:len(read)/2]
					for shingle in getDocShingles(leftpart, k):
						if shingle not in shinglesPos:
							shinglesPos[shingle] = pos
							pos += 1
					rightpart = read[len(read)/2:]
					for shingle in getDocShingles(rightpart, k):
						if shingle not in shinglesPos:
							shinglesPos[shingle] = pos
							pos += 1
					if len(leftpart) < k:
						logprint(log, False,
								 "ERROR: k larger than part length\n",
								 "		 Pick k smaller than", len(leftpart)
						)
						sys.exit()
					read = ""
			# Concatenate multi-line sequences into one string
			else:
				read += line.strip().upper()
		# Compute shingles from last read
		if read != "":
			# Splits the string into two parts
			leftpart = read[:len(read)/2]
			for shingle in getDocShingles(leftpart, k):
				if shingle not in shinglesPos:
					shinglesPos[shingle] = pos
					pos += 1
			rightpart = read[len(read)/2:]
			for shingle in getDocShingles(rightpart, k):
				if shingle not in shinglesPos:
					shinglesPos[shingle] = pos
					#pos += 1

		logprint(log, False, "Finished computation of shingles table in",
				 (time.clock() - tim) / 60, "minutes")
		logprint(log, False, "Number of shingles:", len(shinglesPos))
		logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
		return shinglesPos


def computeShinglesSet(fasta_file, shingles, k, log):
	"""
	Computes the set of all k-shingles (k-mers) in all reads.
	"""
	if not fasta_file:
		return shingles
	logprint(log, True, "Computing set of all shingles...")
	tim = time.clock()
	with open(fasta_file, "rU") as fasta_file:
		#shingles = set()
		read = ""
		for line in fasta_file:
			# If line starts with ">", which indicates end of a sequence
			# Then, append it to list of reads
			if line.startswith(">"):
				if read != "":
					# Splits the string into two parts
					leftpart = read[:len(read)/2]
					for shingle in getDocShingles(leftpart, k):
						shingles.add(shingle)
					rightpart = read[len(read)/2:]
					for shingle in getDocShingles(rightpart, k):
						shingles.add(shingle)
					read = ""
			# Concatenate multi-line sequences into one string
			else:
				read += line.strip().upper()
		# Compute shingles from last read
		if read != "":
			# Splits the string into two parts
			leftpart = read[:len(read)/2]
			for shingle in getDocShingles(leftpart, k):
				shingles.add(shingle)
			rightpart = read[len(read)/2:]
			for shingle in getDocShingles(rightpart, k):
				shingles.add(shingle)

		logprint(log, False, "Finished shingling in", (time.clock() - tim) /
				 60, "minutes")
		logprint(log, False, "Number of shingles:", len(shingles))
		logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
		return shingles


def getAllReads(fasta_file, log, parts=True):
	if fasta_file:
		reads = []
		ext = fasta_file.rsplit(".", 1)[1]
		if ext == "fasta" or ext == "fa":
			with open(fasta_file, "r") as f:
				reads = parseFasta(fasta_file, f, log, parts)
		elif ext == "fastq":
			with open(fasta_file, "r") as f:
				reads = parseFastq(fasta_file, f, log, parts)
		elif ext == "gz" or ext == "tar":
			tar = tarfile.open(fasta_file, "r")
			for file in tar.getmembers():
				ext = file.name.rsplit(".", 1)
				if len(ext) > 1:
					ext = file.name.rsplit(".", 1)[1]
				else:
					continue
				print ext
				if ext == "fasta":
					f = tar.extractfile(file)
					reads += parseFasta(fasta_file, f, log, parts)
				elif ext == "fastq":
					f = tar.extractfile(file)
					reads += parseFastq(fasta_file, f, log, parts)
		else:
			print "File format", "*"+ext+"*", "is NOT supported. Use FASTA",
			print "or FASTQ, possibly compressed as tar.gz."
			sys.exit()

		return reads

	else:
		return []


def parseFastq(filename, fasta_f, log, parts=True):
	read = ""
	tim = time.clock()
	logprint(log, False, "Collecting reads from file", filename)
	reads = []
	seqs = 0
	seq_line = False
	for line in fasta_f:
		# If line starts with ">", which indicates end of a sequence, append it to list of reads
		if line.startswith("@"):
			seq_line = True
		elif line.startswith("+"):
			if read != "":
				seqs += 1
				if parts:
					# Splits the string into two parts
					leftpart = read[:len(read)/2]
					rightpart = read[len(read)/2:]
					reads.append(leftpart)
					reads.append(rightpart)
				else:
					reads.append(read)
				read = ""
				seq_line = False
		elif seq_line:
			# Concatenate multi-line sequences into one string
			read += line.strip().upper()

	logprint(log, False, "Finished reading in", (time.clock() - tim) / 60, "minutes")
	logprint(log, False, "Found", seqs, "sequences in fasta file")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	return reads


def parseFasta(filename, fasta_f, log, parts=True):
	"""
	Extract the reads (DNA sequences) from the given fasta file.
	Splits the reads into two equally (for now) sized parts.
	"""
	read = ""
	tim = time.clock()
	logprint(log, False, "Collecting reads from file", filename)
	reads = []
	seqs = 0
	for line in fasta_f:
		# If line starts with ">", which indicates end of a sequence,
		# append it to list of reads
		if line.startswith(">"):
			if read != "":
				seqs += 1
				if parts:
					# Splits the string into two parts
					leftpart = read[:len(read)/2]
					rightpart = read[len(read)/2:]
					reads.append(leftpart)
					reads.append(rightpart)
				else:
					reads.append(read)
				read = ""
		# Concatenate multi-line sequences into one string
		else:
			read += line.strip().upper()
	if read != "":
		seqs += 1
		if parts:
			leftpart = read[:len(read)/2]
			rightpart = read[len(read)/2:]
			reads.append(leftpart)
			reads.append(rightpart)
		else:
			reads.append(read)

	logprint(log, False, "Finished reading in",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, False, "Found", seqs, "sequences in fasta file")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	return reads


def getPartsFromFile(fasta_file, log):
	"""
	Makes a generator object of all left- and right-parts of all reads
	in the given fasta file.
	"""
	if fasta_file:
		with open(fasta_file, "r") as fasta_file:
			read = ""
			for line in fasta_file:
				# If line starts with ">", which indicates end of a sequence, append it to list of reads
				if line.startswith(">"):
					if read != "":
						# Splits the string into two parts
						leftpart = read[:int(len(read)*leftPartRatio)]
						yield leftpart
						rightpart = read[int(len(read)*rightPartRatio):]
						yield rightpart
						read = ""
				# Concatenate multi-line sequences into one string
				else:
					read += line.strip().upper()
			if read != "":
				leftpart = read[:int(len(read)*leftPartRatio)]
				yield leftpart
				rightpart = read[int(len(read)*rightPartRatio):]
				yield rightpart


def isPrime(n):
	"""
	Checks if the given number (n) is a prime number.
	"""
	if n == 2 or n == 3: return True
	if n < 2 or n % 2 == 0: return False
	if n < 9: return True
	if n % 3 == 0: return False
	r = int(n**0.5)
	f = 5
	while f <= r:
		if n % f == 0: return False
		if n % (f+2) == 0: return False
		f += 6
	return True


def getPrime(offset):
	"""
	Finds the first prime number higher than a given offset.
	"""
	start = random.randrange(100)
	# print "start", start
	offset += start
	if offset % 2 == 0:
		offset += 1
	while True:
		if isPrime(offset):
			return offset
		offset += 2


# ************************************************************************** #
#																			 #
#						  Locality Sensitive Hashing						 #
#																			 #
# ************************************************************************** #
def doWork(tup, b=None, q=None):
	if b > 12:
		#time.sleep(9000)
		#time.sleep(b)
		time.sleep(18000)
	if b != None:
		normal, diseased, shingles, k, rows, min_alg, bands, p = tup
	else:
		normal, diseased, shingles, k, rows, min_alg, b, bands, p = tup
	print b
	seqs = getAllReads(normal, None) + getAllReads(diseased, None)
	#r = redis.StrictRedis()
	buckets = dict()
	#candidatePairs = dict()
	num = minhashing(normal, diseased, shingles, buckets, k, rows,
			   min_alg, b, bands, p, None)
	numPairs = lshBandMultiCore(buckets, b, seqs, None, q)

	#numPairs = 0
	#for key in candidatePairs:
		#r.rpush(key, candidatePairs[key])
		#pre = "redis-cli SADD "+str(key)+" "
		#cmd = pre+" ".join(map(str,candidatePairs[key]))
		#os.system(cmd+" > /dev/null")
		#numPairs += len(candidatePairs[key])
		#val = ','.join(map(str, candidatePairs[key]))
		#r.sadd(key,candidatePairs[key])

	# lst = [i for i in xrange(1000)]
	# for i in xrange(100):
	# 	r.rpush(i, lst)
	print "Num buckets", len(buckets)

	#exportCandidatePairs(candidatePairs, filename, None, num)
	#multiSeqAlign(normal, diseased, b, bands, prefix, suffix, num)

	#return buckets
	#sys.exit()
	return None


def runLSH(normal, diseased, bands, rows, k, seed, minhash_alg, test, log, multiProcessing, pool):
	"""
	Minhash algorithms:
		pre-computed hash functions:
			1. First hash
			2. Through whole matrix (according to the book)
			3. Through all documents shingles
		Ongoing hash functions in the form ((ax+b) mod p) mod N:
			4. First hash
			5. Through whole matrix (according to the book)
			6. Through all documents shingles
	"""
	# Check if files are provided
	if normal or diseased:
		tim = time.clock()
		random.seed(seed)
		candidatePairs = dict()

		# getAllReads(normal, log)
		# sys.exit()

		if minhash_alg == 7:
			p = getPrime(4**k)
			shingles = None
		# Computes table of all k-shingles and their position
		elif minhash_alg == 3 or minhash_alg == 6:
			# shingles = computeShinglesTable(fasta_file, k, log)
			shingles = dict()
			shingles = computeShinglesTable(normal, shingles, k, log)
			shingles = computeShinglesTable(diseased, shingles, k, log)
			# p = getPrime(len(shingles))
			p = 0
		# Computes set of all k-shingles
		else:  # minhash alg 1, 2, 4 or 5
			# shingles = computeShinglesSet(fasta_file, k, log)
			shingles = set()
			shingles = computeShinglesSet(normal, shingles, k, log)
			shingles = computeShinglesSet(diseased, shingles, k, log)
			shingles = list(shingles)
			#p = getPrime(len(shingles))
			p = 0
		# Use Locality-Sensitive Hashing to compute for each bands the buckets
		# with similar documents (reads) obtained by minhashing each read.
		if not multiProcessing:
			# seqs = getAllReads(normal, None) + getAllReads(diseased, None)
			for b in xrange(bands):
				buckets = dict()
				minhashing(normal, diseased, shingles, buckets, k, rows,
						   minhash_alg, b, bands, p, log)
				logprint(log, False, "Number of buckets", len(buckets))
				lshBand(buckets, b, candidatePairs, log)

				# Stop if memory limit reached
				# limit = memory_usage_resource()
				# if int(limit) > 300000:
				# 	break

		if multiProcessing:
			if pool:
				numProcs = 3
				start = numProcs if numProcs <= bands else bands
				prev_start = 0
				stop = bands
				while True:
					params = []
					for b in xrange(prev_start, start):
						params.append( (normal, diseased, shingles, k, rows,
										minhash_alg, b, bands, p) )
					prev_start = start
					results = pool.map(doWork, params)
					tim = time.clock()
					logprint(log, False, "Combining candidate pairs...")
					for tempDict in results:
					#for _ in xrange(bands):
						#tempDict = results.get()
						for key in tempDict:
							if key in candidatePairs:
								for item in tempDict[key]:
									candidatePairs[key].add(item)
							else:
								candidatePairs[key] = set(tempDict[key])
					logprint(log, False,
							 "Finished combining candidate pairs in",
						 	 (time.clock() - tim) / 60, "minutes")
					logprint(log, True, "Memory usage (in mb):",
							 memory_usage_resource())

					start += numProcs
					if start > bands:
						start = bands
					if prev_start == bands:
						break
			else:
				q = Queue()
				#seqs = getAllReads(normal, log) + getAllReads(diseased, log)
				params = (normal, diseased, shingles, k, rows,
						   minhash_alg, bands, p)

				numProcs = 25
				start = numProcs if numProcs <= bands else bands
				prev_start = 0
				stop = bands
				count = 0
				while True: #start < stop+1:

					#processes = []
					for b in xrange(prev_start, start):
						p = Process(target=doWork, args=(params, b, q, ))
						p.start()
						#processes.append(p)
					prev_start = start
					# pipe = r.pipeline()

					tim = time.clock()
					logprint(log, False, "Combining candidate pairs...")
					#d = shelve.open("shelveDBs/cache0")
					results = 0
					while count < start:
						key, mates = q.get()
						#time.sleep(0.1)
						if key == -1:
							count += 1
							print "Done bands:", count

						else:
							if key not in candidatePairs:
								candidatePairs[key] = set()
							for mate in mates:
								candidatePairs[key].add(mate)
							# key = str(key)
							# if d.has_key(key):
							# 	temp_set = d[key]
							# 	for mate in mates:
							# 		temp_set.add(mate)
							# 	d[key] = temp_set
							# else:
							# 	d[key] = set(mates)

					start += numProcs
					if start > bands:
						start = bands
					# for key in candidatePairs:
					# 	print key, candidatePairs[key]

					# shelvesList = [shelve.open("shelveDBs/cache"+str(b), "r") for b
					# 	 		   in xrange(1,bands)]
					# for i in d.keys():
					# 	finalSet = set(d[str(i)])
					# 	#for b in xrange(1, bands):
					# 		#d_temp = shelve.open("shelveDBs/cache"+str(b), "r")
					# 	for d_temp in shelvesList:
					# 		lst_temp = d_temp[i]
					# 		#d_temp.close()
					# 		for item in lst_temp:
					# 			finalSet.add(item)
					# 	d[str(i)] = list(finalSet)


					#d.close()
					logprint(log, False,
							 "Finished combining candidate pairs in",
							 (time.clock() - tim) / 60, "minutes")
					logprint(log, True, "Memory usage (in mb):",
							 memory_usage_resource())

					# q.close()
					# for p in processes:
					# 	p.join()
					if prev_start == bands:
						break


		logprint(log, False, "\nNumber of unique candidate pairs",
				 sum(len(candidatePairs[i]) for i in candidatePairs)/2)
		logprint(log, False, "Finished LSH in",
				 (time.clock() - tim) / 60, "minutes")
		logprint(log, True, "Memory usage (in mb):", memory_usage_resource(),
				 "\n")

		# If benchmarking different k-values
		if test == 3:
			return (time.clock() - tim) / 60, memory_usage_resource(), \
					sum(len(candidatePairs[i]) for i in candidatePairs)/2

		return candidatePairs

	else:
		logprint(log, True, "ERROR: NO FASTA FILES OR IMPORT FILE PROVIDED")
		sys.exit()


def getDocShingles(dna, k, asSet=True):
	"""
	Computes all shingles of size k in a document (dna sequence)
	"""
	if asSet:
		shingles = {dna[i:i+k] for i in xrange(len(dna)-k+1)}
	else: # as list
		shingles = [toBase10(dna[i:i+k]) for i in xrange(len(dna)-k+1)]
	return shingles


def toBase10(seq):
	"""Compute the number from base 4 to base b."""
	#digits = []
	n = 0
	for s in seq:
		if s == "A":
			i = 0
		elif s == "C":
			i = 1
		elif s == "G":
			i = 2
		elif s == "T":
			i = 3
		n = 4 * n + i
	return n


def minhashing(normal, diseased, kmers, buckets, k, rows, minhash_alg, bn, bs, p, log):
	tim = time.clock()
	logprint(log, True, "Minhashing...")
	# random.seed(bn)

	idx = 0
	if minhash_alg < 4:
		hashfuncs = computeHashFunctions(rows, kmers, log)
	else:
		if minhash_alg == 7:
			n = 4**k
		else:
			n = len(kmers)
		p = getPrime(n) # temp
		a = [random.randrange(1, n) for i in xrange(rows)]
		b = [random.randrange(n) for i in xrange(rows)]
	for part in getPartsFromFile(normal, log):
		if minhash_alg == 6:  # Default minhash alg 6
			minhashing_alg6(part, idx, kmers, buckets, k, rows, p, a, b, n)
		elif minhash_alg == 7:
			minhashing_alg7(part, idx, kmers, buckets, k, rows, p, a, b, n)
		elif minhash_alg == 1:
			minhashing_alg1(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 2:
			minhashing_alg2(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 3:
			minhashing_alg3(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 4:
			minhashing_alg4(part, idx, kmers, buckets, k, rows, p, a, b)
		elif minhash_alg == 5:
			minhashing_alg5(part, idx, kmers, buckets, k, rows, p, a, b)

		idx += 1

		if idx % printMinhashProcess == 0:
			logprint(log, True, "Band", bn+1, "of", str(bs)+":",
					 "Processed", idx, "documents in",
					 (time.clock() - tim) / 60, "minutes")

	global secondSample
	if secondSample == 0:
		secondSample = idx

	for part in getPartsFromFile(diseased, log):
		if minhash_alg == 6:  # Default minhash alg 6
			minhashing_alg6(part, idx, kmers, buckets, k, rows, p, a, b, n)
		elif minhash_alg == 7:
			minhashing_alg7(part, idx, kmers, buckets, k, rows, p, a, b, n)
		elif minhash_alg == 1:
			minhashing_alg1(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 2:
			minhashing_alg2(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 3:
			minhashing_alg3(part, idx, kmers, buckets, k, rows, hashfuncs)
		elif minhash_alg == 4:
			minhashing_alg4(part, idx, kmers, buckets, k, rows, p, a, b)
		elif minhash_alg == 5:
			minhashing_alg5(part, idx, kmers, buckets, k, rows, p, a, b)

		idx += 1

		if idx % printMinhashProcess == 0:
			logprint(log, True, "Band", bn+1, "of", str(bs)+":",
					 "Processed", idx, "documents in",
					 (time.clock() - tim) / 60, "minutes")

	logprint(log, False, "Finished minhashing in",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	return idx


def lshBand(buckets, b, candidatePairs, log):
	tim = time.clock()
	logprint(log, True, "Running LSH and finding similar pairs...")
	numPairsUnique = 0
	b += 1
	naiveSim = 0.97
	total = 0
	o = 33  # overlap
	skippedBuckets = 0
	maxBucket = 0
	numBucketsPrintet = 0
	# logprint(log, False, "Bucket sizes:")
	for bucket in buckets:
		# if len(buckets[bucket]) > 1:
		# 	logprint(log, False, str(bucket)+":", len(buckets[bucket]))
		# 	if len(buckets[bucket]) > maxBucket:
		# 		maxBucket = len(buckets[bucket])
		# 	numBucketsPrintet += 1
		if len(buckets[bucket]) > 500:
			skippedBuckets += 1
			continue
		for i in xrange(len(buckets[bucket])):
			id1 = buckets[bucket][i]
			for j in xrange(i+1, len(buckets[bucket])):
				id2 = buckets[bucket][j]
				if id1 % 2 == 0 and id2 % 2 == 1:
					if id1 + 1 != id2:
						# naive = globalAlignment(seqs[id1], seqs[id2], o)
						# if naive >= naiveSim:
						# 	if id1 in candidatePairs:
						# 		candidatePairs[id1].add(id2)
						# 	else:
						# 		candidatePairs[id1] = set([id2])
						# 	if id2 in candidatePairs:
						# 		candidatePairs[id2].add(id1)
						# 	else:
						# 		candidatePairs[id2] = set([id1])
						# 	numPairsUnique += 1
						# total += 1
						if id1 in candidatePairs:
							candidatePairs[id1].add(id2)
						else:
							candidatePairs[id1] = set([id2])
						if id2 in candidatePairs:
							candidatePairs[id2].add(id1)
						else:
							candidatePairs[id2] = set([id1])
						numPairsUnique += 1
						# if id1 in candidatePairs:
						# 	candidatePairs[id1].append(id2)
						# else:
						# 	candidatePairs[id1] = [id2]
						# if id2 in candidatePairs:
						# 	candidatePairs[id2].append(id1)
						# else:
						# 	candidatePairs[id2] = [id1]

				elif id1 % 2 == 1 and id2 % 2 == 0:
					if id1 - 1 != id2:
						# naive = globalAlignment(seqs[id1], seqs[id2], o)
						# if naive >= naiveSim:
						# 	if id1 in candidatePairs:
						# 		candidatePairs[id1].add(id2)
						# 	else:
						# 		candidatePairs[id1] = set([id2])
						# 	if id2 in candidatePairs:
						# 		candidatePairs[id2].add(id1)
						# 	else:
						# 		candidatePairs[id2] = set([id1])
						# 	numPairsUnique += 1
						# total += 1

						if id1 in candidatePairs:
							candidatePairs[id1].add(id2)
						else:
							candidatePairs[id1] = set([id2])
						if id2 in candidatePairs:
							candidatePairs[id2].add(id1)
						else:
							candidatePairs[id2] = set([id1])
						numPairsUnique += 1

						# if id1 in candidatePairs:
						# 	candidatePairs[id1].append(id2)
						# else:
						# 	candidatePairs[id1] = [id2]
						# if id2 in candidatePairs:
						# 	candidatePairs[id2].append(id1)
						# else:
						# 	candidatePairs[id2] = [id1]


	logprint(log, True, "Number of buckets in band", str(b)+":", len(buckets))
	logprint(log, False, "Skipped buckets:", skippedBuckets)
	numPairs = 0
	for bucket in buckets:
		# print buckets[bucket]
		numPairs += len(buckets[bucket]) * (len(buckets[bucket])-1) / 2
	logprint(log, False, "Number of candidate pairs in band", str(b)+":",
			 numPairs)
	logprint(log, True, "Number of unique candidate pairs in band",
			 str(b)+":", numPairsUnique)
	logprint(log, False, "Total number of candidate pairs:", total)
	# logprint(log, True, "Ratio:", float(numPairsUnique) / total)

	# temp
	#logprint(log, False, "Number of buckets:", len(buckets))
	#logprint(log, False, "Max bucket size:", maxBucket)
	#logprint(log, False, "Num buckets sizes printet:", numBucketsPrintet)

	# print "Finished LSH for band", b, "in", (time.clock() - tim) / 60, \
	#		"minutes"

	# print len(candidatePairs)

	return None


def lshBandMultiCore(buckets, b, seqs, log, r=None):
	tim = time.clock()
	logprint(log, True, "Running LSH and finding similar pairs...")
	numPairsUnique = 0
	#d = shelve.open("shelveDBs/cache"+str(b), "c")
	b += 1
	naiveSim = 0.97
	total = 0
	o = 33  # overlap
	skippedBuckets = 0

	# SET VERSION
	if setVersion:
		pipe = r.pipeline()
		for bucket in buckets:
			for i in xrange(len(buckets[bucket])):
				id1 = buckets[bucket][i]
				for j in xrange(i+1, len(buckets[bucket])):
					id2 = buckets[bucket][j]
					if id1 % 2 == 0 and id2 % 2 == 1:
						if id1 + 1 != id2:
							pipe.sadd(id1,id2)
							pipe.sadd(id2,id1)
							numPairsUnique += 1
					elif id1 % 2 == 1 and id2 % 2 == 0:
						if id1 - 1 != id2:
							pipe.sadd(id1,id2)
							pipe.sadd(id2,id1)
							numPairsUnique += 1
			pipe.execute()

	# LIST VERSION
	elif r:
		for bucket in buckets:
			if len(buckets[bucket]) > 500:
				skippedBuckets += 1
				continue
			for i in xrange(len(buckets[bucket])):
				id1 = buckets[bucket][i]
				mates = []
				for j in xrange(len(buckets[bucket])):
					id2 = buckets[bucket][j]
					if id1 % 2 == 0 and id2 % 2 == 1:
						if id1 + 1 != id2:
							naive = globalAlignment(seqs[id1], seqs[id2], o)
							if naive >= naiveSim:
								mates.append(id2)
								numPairsUnique += 1
							total += 1
					elif id1 % 2 == 1 and id2 % 2 == 0:
						if id1 - 1 != id2:
							naive = globalAlignment(seqs[id1], seqs[id2], o)
							if naive >= naiveSim:
								mates.append(id2)
								numPairsUnique += 1
							total += 1
				r.put((id1, mates))
		print "Dropping last"
		r.put((-1, -1))
		print "last dropped"

	else:
		candidatePairs = dict()
		for bucket in buckets:
			if len(buckets[bucket]) > 500:
				skippedBuckets += 1
				continue
			for i in xrange(len(buckets[bucket])):
				id1 = buckets[bucket][i]
				for j in xrange(i+1, len(buckets[bucket])):
					id2 = buckets[bucket][j]
					if id1 % 2 == 0 and id2 % 2 == 1:
						if id1 + 1 != id2:
							naive = globalAlignment(seqs[id1], seqs[id2], o)
							if naive >= naiveSim:
								if id1 in candidatePairs:
									candidatePairs[id1].append(id2)
								else:
									candidatePairs[id1] = [id2]
								if id2 in candidatePairs:
									candidatePairs[id2].append(id1)
								else:
									candidatePairs[id2] = [id1]
								numPairsUnique += 1
							total += 1
					elif id1 % 2 == 1 and id2 % 2 == 0:
						if id1 - 1 != id2:
							naive = globalAlignment(seqs[id1], seqs[id2], o)
							if naive >= naiveSim:
								if id1 in candidatePairs:
									candidatePairs[id1].append(id2)
								else:
									candidatePairs[id1] = [id2]
								if id2 in candidatePairs:
									candidatePairs[id2].append(id1)
								else:
									candidatePairs[id2] = [id1]
								numPairsUnique += 1
							total += 1
		#return candidatePairs

		#pipe = r.pipeline()
		# for bucket in buckets:
		# 	leftParts = []
		# 	rightParts = []
		# 	for item in buckets[bucket]:
		# 		if item % 2 == 0:
		# 			leftParts.append(item)
		# 		else:
		# 			rightParts.append(item)
		# 	numPairsUnique += len(leftParts)*len(rightParts)
		# 	for key in leftParts:
		# 		#pipe.rpush(key, rightParts)
		# 		#print key, rightParts
		# 		#d[str(key)] = rightParts
		# 		r.put((key, rightParts))
		# 	for key in rightParts:
		# 		#pipe.rpush(key, leftParts)
		# 		#print key, leftParts
		# 		#d[str(key)] = leftParts
		# 		r.put((key, leftParts))
			#pipe.execute()

	# for key in d.keys():
	# 	print key, d[key]

	logprint(log, False, "Skipped buckets:", skippedBuckets)
	logprint(log, False, "Number of buckets in band", str(b)+":",len(buckets))
	logprint(log, False, "Number of unique candidate pairs in band",
			 str(b)+":", numPairsUnique)
	logprint(log, False, "Total number of candidate pairs:", total)
	logprint(log, True, "Ratio:", float(numPairsUnique) / total)

	#d.close()

	# print "Finished LSH for band", b, "in", (time.clock() - tim) / 60, \
	#		"minutes"

	# print len(candidatePairs)

	if not setVersion and not r:
		return candidatePairs
	return numPairsUnique


# ************************************************************************** #
#																			 #
#							 Minhashing algorithms							 #
#																			 #
# ************************************************************************** #
def minhashing_alg1(dna, idx, shingles, buckets, k, rows, hashfuncs):
	"""
	Uses pre-computed hash functions and first hash
	"""
	# Create minhash signatures as described in chapter 3 of the book Massive
	# Data Mining
	# Find signature for each document
	docShingles = getDocShingles(dna, k)
	signature = [None for i in xrange(rows)]
	# For each row in the 'character matrix'
	for sigPos in xrange(rows):
		for i, h in enumerate(hashfuncs[sigPos]):
			if shingles[h] in docShingles:
				signature[sigPos] = i
				break

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


def minhashing_alg2(dna, idx, shingles, buckets, k, rows, hashfuncs):
	"""
	Uses pre-computed hashFuncs and runs through each hash (whole matrix)
	and saves smallest value which is in doc
	"""
	# Create minhash signatures as described in chapter 3 of the book Massive
	# Data Mining
	# Find signature for each document
	docShingles = getDocShingles(dna, k)
	signature = [None for i in xrange(rows)]
	# For each row in the 'character matrix'
	for r in xrange(len(shingles)):
		# If the shingle is in the document, then
		if shingles[r] in docShingles:
			# Find the 'first' shingle relative to each permutation
			for i in xrange(rows):
				if signature[i] is None or signature[i] > hashfuncs[i][r]:
					signature[i] = hashfuncs[i][r]

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


def minhashing_alg3(dna, idx, shingles, buckets, k, rows, hashfuncs):
	"""
	Uses pre-computed hashFuncs and table to find original shingle position,
	then find new shingle with smallest position in hash function.
	"""
	# Create minhash signatures as described in chapter 3 of the book Massive
	# Data Mining
	# Find signature for each document
	signature = []
	docShingles = getDocShingles(dna, k)
	numShingles = len(shingles)
	for h in hashfuncs:
		minVal = numShingles+1
		for shingle in docShingles:
			pos = shingles[shingle]
			if h[pos] < minVal:
				minVal = h[pos]
		signature.append(minVal)

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


def minhashing_alg4(dna, idx, shingles, buckets, k, rows, p, a, b):
	"""
	Uses hash functions in the form ((a*pos+b) mod p) mod N,
	where a and b random integers and p is prime and p > N.
	Uses the first hash strategy.
	"""
	# Create minhash signatures as described in chapter 3 of the book Massive
	# Data Mining
	# Find signature for each document
	docShingles = getDocShingles(dna, k)
	signature = []
	numShingles = len(shingles)
	# For each row in the 'character matrix'
	for sigPos in xrange(rows):
		for i in xrange(numShingles):
			val = ((a[sigPos]*i+b[sigPos]) % p) % numShingles
			if shingles[val] in docShingles:
				signature.append(i)
				break
		if len(signature)-1 != sigPos:
			# Ugly fix, only happens when hashfunction creates collisions
			signature.append(random.randrange(numShingles))
			#print signature[-1]

	if len(signature) == rows:
		key = ','.join(map(str, signature))
		if key in buckets:
			buckets[key].append(idx)
		else:
			buckets[key] = [idx]
	else:
		print "FUCK MY LIFE"
		sys.exit()


def minhashing_alg5(dna, idx, shingles, buckets, k, rows, p, a, b):
	"""
	Uses hash functions in the form ((a*pos+b) mod p) mod N,
	where a and b random integers and p is prime and p > N.
	Runs through each hash (whole matrix) and saves the smallest value which
	exist in the dna.
	"""
	# Create minhash signatures as described in chapter 3 of the book Massive
	# Data Mining
	# Find signature for each document
	docShingles = getDocShingles(dna, k)
	signature = [None for i in xrange(rows)]
	numShingles = len(shingles)
	# For each row in the 'character matrix'
	for r in xrange(numShingles):
		# If the shingle is in the document, then
		if shingles[r] in docShingles:
			# Find the 'first' shingle relative to each permutation
			for i in xrange(rows):
				pos = ((a[i]*r+b[i]) % p) % numShingles
				if signature[i] is None or signature[i] > pos:
					signature[i] = pos

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


def minhashing_alg6(dna, idx, shingles, buckets, k, rows, p, a, b, n):
	"""
	DEFAULT MINHASH ALGORITHM
	Uses hash functions in the form ((a*pos+b) mod p) mod N,
	where a and b random integers and p is prime and p > N.
	Computes original position of shingle by finding all shingles and
	enumerating them, then store them in a table for fast look up.
	Table is called shingles.
	"""
	# Find signature for each document
	signature = []
	#signature = array.array('l')
	docShingles = getDocShingles(dna, k)
	#numShingles = len(shingles)
	for i in xrange(rows):
		minVal = n
		for shingle in docShingles:
			pos = shingles[shingle]
			val = ((a[i]*pos+b[i]) % p) % n
			if val < minVal:
				minVal = val
		signature.append(minVal)
	# print signature

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


def minhashing_alg7(dna, idx, shingles, buckets, k, rows, p, a, b, n):
	"""
	Uses hash functions in the form ((a*pos+b) mod p) mod N,
	where a and b random integers and p is prime and p > N.
	Computes original position of shingle by finding all shingles and
	enumerating them, then store them in a table for fast look up.
	Table is called shingles.
	"""
	# Find signature for each document
	signature = []
	#signature = array.array('l')
	docShingles = getDocShingles(dna, k, False)
	for i in xrange(rows):
		minVal = n
		for pos in docShingles:
			#pos = shingles[shingle]
			val = ((a[i]*pos+b[i]) % p) # % n
			if val < minVal:
				minVal = val
		signature.append(minVal)
	#print signature

	key = ','.join(map(str, signature))
	if key in buckets:
		buckets[key].append(idx)
	else:
		buckets[key] = [idx]


# ************************************************************************** #
#																			 #
#							  Similarity checkers							 #
#																			 #
# ************************************************************************** #
def pairsFoundByLSH(normal, diseased, candidatePairs, k, b, r, log):
	"""
	Check which pairs, with similarity above the threshold t, are found by
	LSH and which are missed
	"""
	seqsNormal = getAllReads(normal, log)
	global secondSample
	secondSample = len(seqsNormal)
	seqsDiseased = getAllReads(diseased, log)
	seqs = seqsNormal + seqsDiseased

	filename = "alligator_20K_new5"
	path = "lshPairsVsAllPairs/"
	f1 = open(path+"naive_pairs_all_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')
	f2 = open(path+"jaccard_sets_pairs_all_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')
	f3 = open(path+"jaccard_bags_pairs_all_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')
	f4 = open(path+"naive_pairs_lsh_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')
	f5 = open(path+"jaccard_sets_pairs_lsh_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')
	f6 = open(path+"jaccard_bags_pairs_lsh_"+filename+"_b_"+str(b)+"_r_"+
			  str(r)+"_k_"+str(k)+".txt", 'w')


	count = 0
	numPairs = len(seqs) * (len(seqs)-1)
	#sims = dict()
	truePairs_naive = set()
	truePairs_sets = set()
	truePairs_bags = set()
	sim_threshold = 0.4
	naive_threshold = 0.97
	doPrint = False
	tim = time.clock()

	truePairs_lsh_naive = set()
	truePairs_lsh_sets = set()
	truePairs_lsh_bags = set()
	# Compute similarities for all pairs
	for i in xrange(0,len(seqs),2):
		for j in xrange(1, len(seqs),2):
			if i+1 != j:
				count += 1
				naive = globalAlignment(seqs[i],seqs[j], 33)
				jaccard_sets = jaccardSim(seqs[i], seqs[j], k)
				jaccard_bags = jaccardSim(seqs[i], seqs[j], k, False)
				#sims[(i,j)] = (naive, jaccard_sets, jaccard_bags)
				if naive > 0.8 and jaccard_sets < 0.4:
					logprint(log, False, seqs[i])
					logprint(log, False, seqs[j])
					logprint(log, True, "Naive:", naive, "Jaccard:",
							 jaccard_sets,"\n")
				if naive > naive_threshold:
					truePairs_naive.add((i,j))
					f1.write(str(i)+","+str(j)+" "+str(naive)+"\n")
					# if i in candidatePairs:
					#	  if j in candidatePairs[i]:
					#		  truePairs_lsh_naive.add((i,j))
					#		  f4.write(str(i)+","+str(j)+" "+str(naive)+"\n")
				if jaccard_sets > sim_threshold:
					truePairs_sets.add((i,j))
					f2.write(str(i)+","+str(j)+" "+str(jaccard_sets)+"\n")
					# if i in candidatePairs:
					#	  if j in candidatePairs[i]:
					#		  truePairs_lsh_sets.add((i,j))
					#		  f5.write(str(i)+","+str(j)+" "+
					#				   str(jaccard_sets)+"\n")
				if jaccard_bags > sim_threshold:
					truePairs_bags.add((i,j))
					f3.write(str(i)+","+str(j)+" "+str(jaccard_bags)+"\n")
					# if i in candidatePairs:
					#	  if j in candidatePairs[i]:
					#		  truePairs_lsh_bags.add((i,j))
					#		  f6.write(str(i)+","+str(j)+" "+
					#				   str(jaccard_bags)+"\n")
				if i in candidatePairs:
					if j in candidatePairs[i]:
						truePairs_lsh_naive.add((i,j))
						f4.write(str(i)+","+str(j)+" "+str(naive)+"\n")
						truePairs_lsh_sets.add((i,j))
						f5.write(str(i)+","+str(j)+" "+str(jaccard_sets)+"\n")
						truePairs_lsh_bags.add((i,j))
						f6.write(str(i)+","+str(j)+" "+str(jaccard_bags)+"\n")
				if doPrint:
					print i,j
					print seqs[i], seqs[j]
					print naive
					print jaccard_sets
					print jaccard_bags
				if count % 500000 == 0:
					logprint(log, False, "Processed", format(count, ',d'),
							 "pairs in", (time.clock() - tim) / 60, "minutes")
					logprint(log, True, "Memory usage (in mb):",
							 memory_usage_resource())
	processing_time = (time.clock() - tim) / 60
	c = "{:,}".format(count).replace(",", ".")
	logprint(log, False, "Processed", c, "pairs in", processing_time,
			 "minutes")
	logprint(log, False, "Memory usage (in mb):", memory_usage_resource())
	# logprint(log, False, "Difference jaccard sets vs naive\n",
	#		   truePairs_sets.difference(truePairs_naive))
	# logprint(log, False, "Difference naive vs jaccard sets\n",
	#		   truePairs_naive.difference(truePairs_sets))
	logprint(log, False, "Number of all pairs:", count)

	# Compute similarites for lsh pairs
	# totalPairs = 0
	# truePairs_lsh_naive = set()
	# truePairs_lsh_sets = set()
	# truePairs_lsh_bags = set()
	# for i in candidatePairs:
	#	  for j in candidatePairs[i]:
	#		  if i % 2 == 0:
	#			  totalPairs += 1
	#			  if sims[(i,j)][0] > sim_threshold:
	#				  truePairs_lsh_naive.add((i,j))
	#				  f4.write(str(i)+","+str(j)+" "+str(sims[(i,j)][0])+"\n")
	#			  if sims[(i,j)][1] > sim_threshold:
	#				  truePairs_lsh_sets.add((i,j))
	#				  f5.write(str(i)+","+str(j)+" "+str(sims[(i,j)][1])+"\n")
	#			  if sims[(i,j)][2] > sim_threshold:
	#				  truePairs_lsh_bags.add((i,j))
	#				  f6.write(str(i)+","+str(j)+" "+str(sims[(i,j)][2])+"\n")
	logprint(log, False, "Naive pairs not found by LSH\n",
			 len(truePairs_naive.difference(truePairs_lsh_naive)))
	logprint(log, False, "Jaccard set pairs not found by LSH\n",
			 len(truePairs_sets.difference(truePairs_lsh_sets)))
	logprint(log, False, "Jaccard bag pairs not found by LSH\n",
			 len(truePairs_bags.difference(truePairs_lsh_bags)))
	# logprint(log, False, "Number of lsh pairs:", totalPairs)


def findSimilarPairs(reads, candidatePairs, k, b, r, m, log):
	"""
	Find candidate pairs that has a similarity above the threshold t
	"""
	p = "coord_output/"
	f1 = open(p+"naive_vs_jaccard_standard_NA19240_b"+str(b)+"_r"+
			  str(r)+"_k"+str(k)+"_m"+str(m)+".txt", 'w')
	f2 = open(p+"naive_NA19240_b"+str(b)+"_r"+str(r)+
			  "_k"+str(k)+"_m"+str(m)+".txt", 'w')
	f3 = open(p+"standard_jaccard_NA19240_b"+str(b)+"_r"+str(r)+
			  "_k"+str(k)+"_m"+str(m)+".txt", 'w')
	f4 = open(p+"pairs_NA19240_b"+str(b)+"_r"+str(r)+
			  "_k"+str(k)+"_m"+str(m)+".txt", 'w')
	f5 = open(p+"naive_and_jaccard_info_NA19240_b"+str(b)+"_r"+
			  str(r)+"_k"+str(k)+"_m"+str(m)+".txt", 'w')

	counter = 0
	logprint(log, True, "Finding similar pairs")
	timer = time.clock()
	numReads = len(reads)
	maxNumPairs = len(candidatePairs) * (len(candidatePairs)-1)
	logprint(log, True, "Number of reads", numReads)
	for doc1 in candidatePairs:
		for doc2 in candidatePairs[doc1]:
			counter += 1
			dna1 = reads[doc1]
			dna2 = reads[doc2]
			if doc1 % 2 == 0:
				bestMatch1 = globalAlignment(dna1, dna2, False)
				bestMatch2 = jaccardSim(dna1, dna2, k)
			else:
				bestMatch1 = globalAlignment(dna2, dna1, False)
				bestMatch2 = jaccardSim(dna2, dna1, k)

			if extraInfo:
				f1.write(str(bestMatch1[0]) + " " + str(bestMatch2) + "\n")
				f2.write(str(bestMatch1[0]) + " " + str(bestMatch1[1]) + " " +
						 str(bestMatch1[2]) + "\n")
				f5.write(str(bestMatch1[0]) + " " + str(bestMatch1[1]) + " " +
						 str(bestMatch1[2]) + " " + str(bestMatch2) + "\n")
			else:
				f1.write(str(bestMatch1) + " " + str(bestMatch2) + "\n")
				f2.write(str(bestMatch1) + "\n")
			f3.write(str(bestMatch2) + "\n")
			f4.write(str(doc1) + " " + str(doc2) + "\n")

			if counter % 500000 == 0:
				logprint(log, False, "Processed", format(counter, ',d'),
						 "pairs in", (time.clock() - timer) / 60, "minutes")
				logprint(log, True, "Memory usage (in mb):",
						 memory_usage_resource())

	processing_time = (time.clock() - timer) / 60
	c = "{:,}".format(counter).replace(",", ".")
	logprint(log, False, "Processed", c, "pairs in", processing_time,
			 "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def globalAlignment(dna1, dna2, t, extraInfo=False):
	"""
	Aligning sequences by using a sliding window approach.
	Returns the best score (matches / seqlength) between the two sequences.
	doc1 is a leftpart and doc2 is a rightpart
	"""
	start = 0
	start2 = 0
	if len(dna1) > len(dna2):
		#start2 = len(dna1) - len(dna2)
		readLen = len(dna2)
	else:
		readLen = len(dna1)
	if extraInfo:
		bestScore = (0, 0, 0)
	else:
		bestScore = 0
	seqLength = readLen-start
	while seqLength > t:
		# print seqLength, bestScore[1]
		matches = 0
		matches2 = 0
		for i in xrange(seqLength):
			# print len(doc1.dna)-start
			if dna1[i] == dna2[i+start]:
				matches += 1
			if dna1[i+start] == dna2[i]:
				matches2 += 1
		# print bestScore
		score = matches / float(seqLength)
		score2 = matches2 / float(seqLength)
		if extraInfo:
			if score > bestScore[0]:
				# print score, bestScore[0]
				# print seqLength, matches, bestScore[1]
				bestScore = (score, matches, seqLength)
				if bestScore[0] == 1.0:
					return bestScore
			if score2 > bestScore[0]:
				bestScore = (score2, matches2, seqLength)
				if bestScore[0] == 1.0:
					return bestScore
		else:
			if score > bestScore:
				bestScore = score
				if bestScore == 1.0:
					return bestScore
			if score2 > bestScore:
				bestScore = score2
				if bestScore == 1.0:
					return bestScore
		start += 1
		seqLength = readLen-start
	return bestScore


def jaccardSim(doc1, doc2, k, jaccard_sets=True):
	"""
	Computing the jaccard similarity.
	Option to use jaccard bag similarity or standard jaccard similarity.
	"""
	# ## standard jaccard sim ## #
	if jaccard_sets:
		shingles1 = getDocShingles(doc1, k)
		shingles2 = getDocShingles(doc2, k)
		intersection = shingles1.intersection(shingles2)
		if len(intersection) == 0:
			return 0
		union = shingles1.union(shingles2)
		return float(len(intersection)) / len(union)

	# ## Bag jaccard sim ## #
	else:
		shingles1 = getDocShingles(doc1, k, False)
		shingles2 = getDocShingles(doc2, k, False)
		counterA = Counter(shingles1)
		counterB = Counter(shingles2)
		intersection = sum((counterA & counterB).values())
		if intersection == 0:
			return 0
		# Definition 1
		#union = len(shingles1) + len(shingles2)# - intersection
		# Definition 2
		union = sum((counterA | counterB).values())
		return float(intersection) / union


def makeSPlot(fasta_file, candidatePairs, k, b, r, alg, log):
	"""
	Method for used for obtaining information about the validity of the pairs
	found by LSH
	"""
	### Used for saving each individual pairs similarity
	# f1 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
	#			"_m"+str(alg)+"_all.txt", 'w')
	# f2 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
	#			"_m"+str(alg)+"_lsh.txt", 'w')
	# f3 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
	#			"_m"+str(alg)+"_rest.txt", 'w')

	reads = getAllReads(fasta_file, log)

	# Create a mapper to look up position in array given a similarity
	n1 = len(reads[0])
	n2 = len(reads[1])
	maxNumShingles = (n1-k+1)+(n2-k+1)
	possibleSims = set([0])
	for i in xrange(1, n1):
		for j in xrange(i, maxNumShingles-i+1):
			possibleSims.add(float(i)/j)
	posSimilarities = dict()
	possibleSims = sorted(list(possibleSims))
	for idx, sim in enumerate(possibleSims):
		posSimilarities[sim] = idx
	# print posSimilarities
	# print len(posSimilarities)

	allPairs = [0 for i in xrange(len(possibleSims))]
	lshPairs = [0 for i in xrange(len(possibleSims))]

	numReads = len(reads)
	pairNum = 0
	timer = time.clock()
	count = 0
	count2 = 0
	process = 0
	# for doc1 in xrange(numReads):
	#	  start = 1
	#	  if doc1 % 2 == 0:
	#		  start = 3
	#	  for doc2 in xrange(doc1+start, numReads, 2):
	for doc1 in xrange(0, numReads, 2):
		for doc2 in xrange(1, numReads, 2):
			if i+1 != j:
				dna1 = reads[doc1]
				dna2 = reads[doc2]
				jaccard = jaccardSim(dna1, dna2, k)

				### Saves the similarity for each pair, and for pairs found by LSH
				# if jaccard > 0:
				#	  f1.write(str(pairNum) + " " + str(jaccard) + "\n")
				#	  rest = True
				#	  if doc1 in candidatePairs:
				#		  if doc2 in candidatePairs[doc1]:
				#			  f2.write(str(pairNum) + " " + str(jaccard) + "\n")
				#			  rest = False
				#			  count += 1
				#	  if rest:
				#		  f3.write(str(pairNum) + " " + str(jaccard) + "\n")
				#		  count2 += 1
				#	  pairNum += 1

				if jaccard > 0:
					pos = posSimilarities[jaccard]
					if doc1 in candidatePairs:
						if doc2 in candidatePairs[doc1]:
							lshPairs[pos] += 1
							count += 1
					allPairs[pos] += 1
					pairNum += 1
				else:
					allPairs[0] += 1


				process += 1
				if process % 500000 == 0:
					logprint(log, True, "Processed", process,
							 "pairs in time:", (time.clock() - timer),
							 "Found", pairNum, "cand. pairs")

	p = "s_plot/"
	f = open(p+"s_shape_info_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
			 "_m"+str(alg)+"_readsfa.txt", 'w')
	for i in xrange(len(allPairs)):
		if allPairs[i] == 0:
			f.write(str(0)+" "+str(0)+" "+str(0)+
					" "+str(possibleSims[i])+"\n")
		else:
			f.write(str(lshPairs[i])+" "+str(allPairs[i])+" "+
				str(float(lshPairs[i]) / allPairs[i])+" "+
				str(possibleSims[i]) + "\n")

	logprint(log, False, "Candidate pairs found by LSH:", count)
	logprint(log, False, "Number of pairs not found by LSH:", pairNum-count)
	logprint(log, True, "Total number of pairs:", pairNum)


def euclideanDistance(doc1, doc2):
	sig1, sig2 = doc1.signature, doc2.signature
	# lol = sum([(x-y)**2 for x, y in itertools.izip(sig1, sig2)])
	# if counter2 < 100:
	#	  print doc1.id, doc2.id
	#	  print sig1
	#	  print sig2
	#	  print lol
	# else:
	#	  sys.exit(0)

	intersum = 0
	for x, y in itertools.izip(sig1, sig2):
		intersum


def longest_common_substring(s1, s2):
	m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
	longest, x_longest = 0, 0
	for x in xrange(1, 1 + len(s1)):
		for y in xrange(1, 1 + len(s2)):
			if s1[x - 1] == s2[y - 1]:
				m[x][y] = m[x - 1][y - 1] + 1
				if m[x][y] > longest:
					longest = m[x][y]
					x_longest = x
			else:
				m[x][y] = 0
	return s1[x_longest - longest: x_longest]


def testLCS(doc1, doc2):
	"""
	Test longest_common_substring method
	"""
	sig1 = ''.join(map(str, doc1.signature))
	sig2 = ''.join(map(str, doc2.signature))
	seq = longest_common_substring(sig1, sig2)
	# print (doc1.id, doc2.id),
	# if counter2 < 100:
	#	  print doc1.id, doc2.id
	#	  print sig1
	#	  print sig2
	#	  print seq
	#	  print
	# else:
	#	  sys.exit(0)
	return seq


def NeedlemanWunsch(doc1, doc2):
	"""
	Sequence alignment using the Needleman-Wunsch algorithm
	"""
	# scores
	match = 1
	mismatch = -1
	indel = -2

	# seqA = "abcd"
	# seqB = "bcd"
	seqA = doc1.dna
	seqB = doc2.dna
	# scoreMatrix = [[0 for x in xrange(len(seqB)+1)] for x in
	#				xrange((len(seqA)+1))]
	scoreMatrix = [[0] * (1 + len(seqB)) for i in xrange(1 + len(seqA))]
	for i in xrange(len(seqA)+1):
		scoreMatrix[i][0] = i*indel

	for j in xrange(len(seqB)+1):
		scoreMatrix[0][j] = j*indel

	for i in xrange(1, len(seqA)+1):
		for j in xrange(1, len(seqB)+1):
			if seqA[i-1] == seqB[j-1]:
				score = scoreMatrix[i-1][j-1] + match
			else:
				score = scoreMatrix[i-1][j-1] + mismatch
			opt1 = scoreMatrix[i-1][j] + indel
			opt2 = scoreMatrix[i][j-1] + indel
			maxi = opt1
			if opt2 > maxi:
				maxi = opt2
			if score > maxi:
				maxi = score
			scoreMatrix[i][j] = maxi
	# print seqA
	# print seqB
	# print scoreMatrix[len(seqA)][len(seqB)]
	# for row in scoreMatrix:
	#	  print row

	# sys.exit(0)


# ************************************************************************** #
#																			 #
#							   Sequence Alignment							 #
#																			 #
# ************************************************************************** #
class AlignedGroup(object):
	"""

	"""
	consensus = []
	#preConsensus = []
	readROffset = 0
	consensusMain = 0  # read_R
	# reads in group coming from left part - Normal sample
	leftPartsN = dict()
	# reads in group coming from left part - Diseased sample
	leftPartsD = dict()
	leftReadsOffset = 0
	maxLeftReadsOffset = 0
	# rightParts = dict()  # reads in group coming from right part
	# reads in group coming from right part and touching pre-consensus
	# preRightParts = dict()
	checkedRightParts = set()
	rightPartGroups = [] # List of RightPart objects
	mismatches = set()

	# Initializer
	def __init__(self, read_R, readROffset, leftReadsOffset):
		self.consensus = []
		#self.preConsensus = []
		self.readROffset = readROffset
		self.consensusMain = read_R
		self.leftPartsN = dict()
		self.leftPartsD = dict()
		self.leftReadsOffset = leftReadsOffset
		self.maxLeftReadsOffset = leftReadsOffset
		#self.rightParts = dict()
		#self.preRightParts = dict()
		self.checkedRightParts = set()
		self.rightPartGroups = []
		self.mismatches = set()


class RightPartGroup(object):
	"""

	"""
	consensus = []
	preConsensus = []
	rightPartsN = dict()
	rightPartsD = dict()
	mismatches = set()

	# Initializer
	def __init__(self, consensus):
		self.consensus = copy.deepcopy(consensus)
		self.preConsensus = []
		self.rightPartsN = dict()
		self.rightPartsD = dict()
		self.mismatches = set()


def print_compare(read_R, misPos, mut, consens, log):
	compareString = ""
	for i in xrange(len(read_R)):
		if i == misPos:
			compareString += mut
		else:
			compareString += read_R[i]
	logprint(log, False, "searchString:", compareString)


def print_fullConsensus(preconsensus, consensus, log=None):
	alphabetSize = 4
	for i in xrange(alphabetSize):
		consensusString = ""
		# pre-consensus
		for j in xrange(len(preconsensus)):
			consensusString += " "
		# for j in xrange(len(preconsensus)):
		# 	if i < len(preconsensus[j]) and \
		# 			preconsensus[j][preconsensus[j].keys()[i]] >= \
		# 			requiredOverlaps:
		# 		consensusString += preconsensus[j].keys()[i]
		# 	else:
		# 		consensusString += " "
		# Space between pre-consensus and consensus
		# if consensusString != "":
		#	  consensusString += " "
		# Get main consensus
		for j in xrange(len(consensus)):
			if i < len(consensus[j]) and \
					consensus[j][consensus[j].keys()[i]] >= requiredOverlaps:
				consensusString += consensus[j].keys()[i]
			else:
				consensusString += " "
		if consensusString.strip() != "":
			print "", consensusString
			if log:
				log.write(" "+consensusString+"\n")
		else:
			continue


def print_alignedGroup(group, rightPartGroup, read_R, seqs, log):
	logprint(log, False, "\nread_R:", read_R)
	logprint(log, False, "Consensus:")
	print_fullConsensus(rightPartGroup.preConsensus,
						rightPartGroup.consensus, log)
	lenPre = len(rightPartGroup.preConsensus)
	newOffset = lenPre + group.readROffset
	# Print normal reads
	if read_R < secondSample:
		logprint(log, False, " " * lenPre,
			seqs[read_R-1]+""+seqs[read_R]+"*")
	for read_L in group.leftPartsN:
		for offset in group.leftPartsN[read_L]:
			logprint(log, False, " " * (newOffset +
				 offset), seqs[read_L]+""+seqs[read_L+1])
	for read_R2 in rightPartGroup.rightPartsN:
		for offset in rightPartGroup.rightPartsN[read_R2]:
			logprint(log, False, " " * (offset + lenPre),
					 seqs[read_R2-1]+""+seqs[read_R2])
	# Print diseased reads
	logprint(log, False, "")
	if read_R >= secondSample:
		logprint(log, False, " " * lenPre,
			seqs[read_R-1]+""+seqs[read_R]+"*")
	for read_L in group.leftPartsD:
		for offset in group.leftPartsD[read_L]:
			logprint(log, False, " " * (newOffset +
				 offset), seqs[read_L]+""+seqs[read_L+1])
	for read_R2 in rightPartGroup.rightPartsD:
		for offset in rightPartGroup.rightPartsD[read_R2]:
			logprint(log, False, " " * (offset + lenPre),
					 seqs[read_R2-1]+""+seqs[read_R2])
	leftParts = group.leftPartsN.keys() + \
				group.leftPartsD.keys()
	logprint(log, False, "Left parts:",
			 sorted(list(leftParts)))
	logprint(log, False, "Number of left parts:",
			 len(leftParts))
	rightParts = rightPartGroup.rightPartsN.keys() + \
				 rightPartGroup.rightPartsD.keys()
	logprint(log, False, "Right parts:",
			 sorted(rightParts))
	logprint(log, False, "Number of right parts:",
			 len(rightParts))
	logprint(log, True, "mismatches:",
			 list(rightPartGroup.mismatches))


def print_leftGroup(group, read_R, seqs, log):
	logprint(log, False, "\nread_R:", read_R)
	logprint(log, False, "Consensus:")
	print_fullConsensus([], group.consensus, log)
	logprint(log, False, "", seqs[read_R-1]+""+seqs[read_R]+"*")
	for read_L in group.leftPartsN:
		for offset in group.leftPartsN[read_L]:
			logprint(log, False, " " * (offset + group.readROffset),
					 seqs[read_L]+""+seqs[read_L+1])
	logprint(log, False, "")
	for read_L in group.leftPartsD:
		for offset in group.leftPartsD[read_L]:
			logprint(log, False, " " * (offset + group.readROffset),
					 seqs[read_L]+""+seqs[read_L+1])


def print_alignedGroups(groups, read_R, seqs, log):
	for group in groups:
		if len(group.rightPartGroups) > 0:
			for rightPartGroup in group.rightPartGroups:
				if len(rightPartGroup.mismatches) > 0:
					logprint(log, False, "\nread_R:", read_R)
					logprint(log, False, "Consensus:")
					# logprint(log, False, "", ''.join(consensus.keys()[0]
					#		  for consensus in group.preConsensus) + #" " +
					#		  ''.join(consensus.keys()[0] for consensus
					#		  in group.consensus))
					print_fullConsensus(rightPartGroup.preConsensus,
										rightPartGroup.consensus, log)
					lenPre = len(rightPartGroup.preConsensus)
					newOffset = lenPre + group.readROffset
					# Print normal reads
					if read_R < secondSample:
						logprint(log, False, " " * lenPre,
							seqs[read_R-1]+""+seqs[read_R]+"*")
					for read_L in group.leftPartsN:
						for offset in group.leftPartsN[read_L]:
							logprint(log, False, " " * (newOffset +
								 offset), seqs[read_L]+""+seqs[read_L+1])
					for read_R2 in rightPartGroup.rightPartsN:
						for offset in rightPartGroup.rightPartsN[read_R2]:
							logprint(log, False, " " * (offset + lenPre),
									 seqs[read_R2-1]+""+seqs[read_R2])
					# Print diseased reads
					logprint(log, False, "")
					if read_R >= secondSample:
						logprint(log, False, " " * lenPre,
							seqs[read_R-1]+""+seqs[read_R]+"*")
					for read_L in group.leftPartsD:
						for offset in group.leftPartsD[read_L]:
							logprint(log, False, " " * (newOffset +
								 offset), seqs[read_L]+""+seqs[read_L+1])
					for read_R2 in rightPartGroup.rightPartsD:
						for offset in rightPartGroup.rightPartsD[read_R2]:
							logprint(log, False, " " * (offset + lenPre),
									 seqs[read_R2-1]+""+seqs[read_R2])
					leftParts = group.leftPartsN.keys() + \
								group.leftPartsD.keys()
					logprint(log, False, "Left parts:",
							 sorted(list(leftParts)))
					logprint(log, False, "Number of left parts:",
							 len(leftParts))
					rightParts = rightPartGroup.rightPartsN.keys() + \
								 rightPartGroup.rightPartsD.keys()
					logprint(log, False, "Right parts:",
							 sorted(rightParts))
					logprint(log, False, "Number of right parts:",
							 len(rightParts))
					logprint(log, False, "mismatches:",
							 list(rightPartGroup.mismatches))
					# logprint(log, False, "Number of mismatches:",
					#		   len(rightPartGroup.mismatches))
					# if len(group.mismatches) > 0:
					#	  sys.exit()
		# else:
		elif len(group.mismatches) > 0:
			logprint(log, False, "\nread_R:", read_R)
			logprint(log, False, "Consensus:")
			# logprint(log, False, "", ''.join(consensus.keys()[0]
			#		   for consensus in group.consensus))
			print_fullConsensus([], group.consensus, log)
			if read_R < secondSample:
				logprint(log, False, "",
						 seqs[read_R-1]+""+seqs[read_R]+"*")
			for read_L in group.leftPartsN:
				for offset in group.leftPartsN[read_L]:
					logprint(log, False, " " *(offset +
							 group.readROffset),
							 seqs[read_L]+""+seqs[read_L+1])
			if read_R >= secondSample:
				logprint(log, False, "",
						 seqs[read_R-1]+""+seqs[read_R]+"*")
			for read_L in group.leftPartsD:
				for offset in group.leftPartsD[read_L]:
					logprint(log, False, " " *(offset +
							 group.readROffset),
							 seqs[read_L]+""+seqs[read_L+1])
			logprint(log, False, "mismatches:", group.mismatches)
			leftParts = group.leftPartsN.keys() + group.leftPartsD.keys()
			logprint(log, False, "Left parts:", sorted(leftParts))
			logprint(log, True, "Number of left parts:", len(leftParts))


def sequenceAlignment(candidatePairs, normal, diseased, log):
	seqsNormal = getAllReads(normal, log)
	global secondSample
	secondSample = len(seqsNormal)
	seqsDiseased = getAllReads(diseased, log)
	seqs = seqsNormal + seqsDiseased

	numAlignedGroups = []
	numRightPartGroups = []

	numMutations1 = 0
	numMutations2 = 0
	numParts = len(candidatePairs) / 4 + 2
	prog = 0
	tim = time.clock()
	for read_R in candidatePairs:
	#for read_R in xrange(1, len(seqs)+1, 20):
		#print read_R
		if read_R < secondSample and read_R % 2 == 1:
		# if read_R % 2 == 1:
		# if read_R == 42535:
		# if read_R == 19:
			if len(candidatePairs[read_R]) > maxCandMates:
				continue

			alignedGroups = []

			# Align left parts
			alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log)

			# Align right parts
			alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log)

			# Analyze the aligned group to find mutations
			if MUTFIND == 1:
				numMutations1 += oldFindMutation(read_R, seqs,
												 alignedGroups,log)
			else:
				numMutations2 += newFindMutation(read_R, seqs,
												 alignedGroups,log)

			# Statistics on number of created groups
			numAlignedGroups.append(len(alignedGroups))
			for group in alignedGroups:
				numRightPartGroups.append(len(group.rightPartGroups))


			#print_alignedGroups(alignedGroups, read_R, seqs, log)
			# sys.exit()
			prog += 1
			if prog % 500 == 0:
				logprint(log, False, "Processed", prog, "of", numParts,
						 "anchor points in", (time.clock()-tim)/60, "minutes")
				logprint(log, True, "Memory usage (in mb):",
						 memory_usage_resource())
				global c1
				logprint(log, False, "left parts aligned:", c1)
				#c1 = 0
				global c2
				logprint(log, True, "right parts aligned:", c2)
				logprint(log, True, "num useful groups:", numMutations1)
				logprint(log, True, "num useful groups:", numMutations2)
				global c3
				logprint(log, False, "positions compared:", c3)
				#c2 = 0

	logprint(log, False, "Finished sequence alignment",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	logprint(log, False, "c1:", c1)
	logprint(log, False, "c2:", c2)
	logprint(log, False, "c3:", c3)
	logprint(log, False, "c4:", c4)
	logprint(log, False, "c5:", c5)
	logprint(log, False, "c6:", c6)
	logprint(log, False, "c7:", c7)
	logprint(log, False, "Too small groups:", c8)
	logprint(log, False, "numReadL:", numreadL)
	logprint(log, False, "numReadR:", numreadR)
	logprint(log, False, "numMutations1:", numMutations1)
	logprint(log, False, "numMutations2:", numMutations2)
	logprint(log, False, "counterAlignGroups:\n", Counter(numAlignedGroups))
	logprint(log, False, "counterRightGroups:\n", Counter(numRightPartGroups))


def multiSeqAlign(seqs, p, pool_size, num, candidatePairs, log):
	#seqs, p, pool_size, num, log = tup
	global p_id
	p_id = p
	log = None
	numMutations1 = 0
	numMutations2 = 0
	numParts = ((num+1) / 2) / pool_size
	prog = 0
	tim = time.clock()
	#r = redis.StrictRedis()
	#r = shelve.open("shelveDBs/cache0", "r")
	for read_R in xrange(p*2+1, secondSample+1, pool_size*2):
		#print read_R

		alignedGroups = []

		# Align left parts
		#alignLeftParts(read_R, seqs, alignedGroups, None, log, r)
		alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log)

		# Align right parts
		#alignRightParts(read_R, seqs, alignedGroups, None, log, r)
		alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log)

		# Analyze the aligned group to find mutations
		if MUTFIND == 1:
			numMutations1 += oldFindMutation(read_R, seqs,
											 alignedGroups, log)
		else:
			numMutations2 += newFindMutation(read_R, seqs,
											 alignedGroups, log)


		#print_alignedGroups(alignedGroups, read_R, seqs, log)
		# sys.exit()
		prog += 1
		if prog % 500 == 0:
			logprint(log, False, "Processed", prog, "of", numParts,
					 "anchor points in", (time.clock()-tim)/60, "minutes")
			logprint(log, True, "Memory usage (in mb):",
					 memory_usage_resource())
			global c1
			logprint(log, False, "left parts aligned:", c1)
			#c1 = 0
			global c2
			logprint(log, True, "right parts aligned:", c2)
			logprint(log, True, "num useful groups:", numMutations1)
			logprint(log, True, "num useful groups:", numMutations2)
			global c3
			logprint(log, False, "positions compared:", c3)
			#c2 = 0

	logprint(log, False, "Finished sequence alignment",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	logprint(log, False, "c1:", c1)
	logprint(log, False, "c2:", c2)
	logprint(log, False, "c3:", c3)
	logprint(log, False, "c4:", c4)
	logprint(log, False, "c5:", c5)
	logprint(log, False, "c6:", c6)
	logprint(log, False, "c7:", c7)
	logprint(log, False, "Too small groups:", c8)
	logprint(log, False, "numReadL:", numreadL)
	logprint(log, False, "numReadR:", numreadR)
	logprint(log, False, "numMutations1:", numMutations1)
	logprint(log, False, "numMutations2:", numMutations2)
	#logprint(log, False, "counterAlignGroups:\n", Counter(numAlignedGroups))
	#logprint(log, False, "counterRightGroups:\n", Counter(numRightPartGroups))
	return p


def initMultiSeqAlign(normal, diseased, candPairs, pool, pool_size, log=None):
	seqsNormal = getAllReads(normal, log)
	global secondSample
	secondSample = len(seqsNormal)
	seqsDiseased = getAllReads(diseased, log)
	seqs = seqsNormal + seqsDiseased
	num = len(seqs)

	#params = [(seqs, p, pool_size, num, log) for p in range(0, pool_size)]
	#params = (seqs, pool_size, num, log)
	#results = pool.map(multiSeqAlign, params)
	for b in xrange(pool_size):
		p=Process(target=multiSeqAlign, args=(seqs, b, pool_size, num,
			 								  candPairs, log, ))
		p.start()
		p.join()
	#multiSeqAlign(params)


def alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log, r=None):
	readROffset = len(seqs[read_R-1])
	checkedLeftParts = set()
	parts = candidatePairs[read_R]
	for read_L in parts:
		if read_L in checkedLeftParts:
			continue
		if read_L < secondSample:
			m = M1  # set to 0 or M1
		else:
			m = M1
		for alignInfo in findAlignment(read_R, read_L, seqs,
										 readROffset, m, log):
			#offset += len(seqs[read_R-1])
			offset, mis, lenCompared = alignInfo
			global numreadL
			numreadL += 1
			#print offset
			newGroup = True
			for group in alignedGroups:
				if fitsInGroup(group, seqs, read_R, read_L, alignInfo,
							   offset, M2):
					# Add read_L to group
					global c1
					c1 += 1
					if read_L < secondSample:
						if read_L in group.leftPartsN:
							group.leftPartsN[read_L].append(offset)
						else:
							group.leftPartsN[read_L] = [offset]
					else:
						if read_L in group.leftPartsD:
							group.leftPartsD[read_L].append(offset)
						else:
							group.leftPartsD[read_L] = [offset]

					# Fits in group, so don't create new group
					newGroup = False
					break
			# If read_L doesn't fit in any group, then create a new
			if newGroup:
				global c7
				c7 += 1
				group = AlignedGroup(read_R, readROffset, offset+readROffset)
				# Start of extra part - offsetExtraPart
				start = len(seqs[read_R]) - offset
				#group.consensusRight = seqs[read_L][start:] + seqs[read_L+1]

				# Add read_L to new the group
				if read_L < secondSample:
					group.leftPartsN[read_L] = [offset]
				else:
					group.leftPartsD[read_L] = [offset]
				group.leftReadsOffset = offset+group.readROffset
				group.mismatches = mis

				# Add anchor point to consensus
				for bp in ''.join((seqs[read_R-1], seqs[read_R])):
					group.consensus.append({bp:1})

				# Add overlapping part of read_L to consensus
				seq_read_L = ''.join((seqs[read_L], seqs[read_L+1]))
				for index in xrange(start):
					i = index + group.readROffset + offset
					bp = seq_read_L[index]
					group.consensus[i][bp] = group.consensus[i].get(bp, 0) + 1

				# print " "*offset, (seqs[read_L]+seqs[read_L+1])
				# print " "*-offset, seqs[read_R]
				# print offset

				# Add the rest of read_L to consensus
				for bp in seq_read_L[start:]:
					group.consensus.append({bp:1})

				# Append new group to the other groups
				alignedGroups.append(group)

			checkedLeftParts.add(read_L)
	# TESTING - Second pass of left-parts
	# for i in xrange(len(alignedGroups)):
	#	  for read_L in alignedGroups[i].leftPartsN:
	#		  for j in xrange(i+1, len(alignedGroups)):
	#			  for offset in alignedGroups[i].leftPartsN[read_L]:
	#				  if fitsInGroup(alignedGroups[j], seqs, read_R, read_L,
	#							   None, offset - group.readROffset, M2):
	#					  if read_L in group.leftPartsN:
	#						  group.leftPartsN[read_L].append(offset)
	#					  else:
	#						  group.leftPartsN[read_L] = [offset]
	#	  for read_L in alignedGroups[i].leftPartsD:
	#		  for j in xrange(i+1, len(alignedGroups)):
	#			  for offset in alignedGroups[i].leftPartsD[read_L]:
	#				  if fitsInGroup(alignedGroups[j], seqs, read_R, read_L,
	#							   None, offset - group.readROffset, M2):
	#					  if read_L in group.leftPartsD:
	#						  group.leftPartsD[read_L].append(offset)
	#					  else:
	#						  group.leftPartsD[read_L] = [offset]


def findAlignment(r_R, r_L, seqs, readROffset, m, log):
	read_R = seqs[r_R]
	read_L = seqs[r_L]
	doPrint = False
	offset = 0
	alignments = 0
	if len(read_R) > len(read_L):
		offset = len(read_R) - len(read_L)
		lengthToCompare = len(read_L)# - overlap*2
	else:
		lengthToCompare = len(read_R)# - overlap*2
	overlapLength = lengthToCompare - overlap*2
	# check for alignment by shifting read_L along read_R
	while lengthToCompare > overlapLength:
		mismatches = set()
		# if log == 2119:
		#	  print "", read_R
		#	  print " "*offset, read_L
		for i in xrange(lengthToCompare):
			if read_R[i+offset] != read_L[i]:
				global c3
				c3 += 1
				mismatches.add(i+offset+readROffset)
				if len(mismatches) > m:
					break
		if len(mismatches) <= m:
			if doPrint:
				print "offset:", offset
				print "", read_R
				print " "*(offset-1), read_L
				print mismatches
			yield offset, mismatches, lengthToCompare
			alignments += 1
			if maxAlignments == alignments:
				break
		offset += 1
		lengthToCompare -= 1


def fitsInGroup(group, seqs, read_R, read_L, alignInfo, offset, m2):
	global c6
	c6 += 1
	lread_R = seqs[read_R]
	lread_L = seqs[read_L]+seqs[read_L+1]
	mismatches = 0
	if alignInfo:
		offset2, mis, lenCompared = alignInfo
	else:
		lenCompared = len(lread_R) - offset
		mis = set()
		for i in xrange(lenCompared):
			if lread_R[i+offset] != lread_L[i]:
				mis.add(i+offset+group.readROffset)
	offset += group.readROffset

	mismatches = set()
	seq_read_R = seqs[read_R-1]+seqs[read_R]
	# lenToCompare = min(len(group.consensus)-offset, len(lread_L))
	# for i in xrange(lenToCompare):
	#	  if len(group.consensus[i+offset]) > 1 or \
	#		 group.consensus[i+offset].keys()[0] != lread_L[i]:
	#		  mismatches.add(i+offset)

	lenToCompare = len(group.consensus)-len(seq_read_R)
	readLLength = len(lread_L) - lenCompared
	# if offset < 50:
	#	  print
	#	  print_fullConsensus([], group.consensus)
	#	  print lenToCompare
	#	  print readLLength
	#	  print lenCompared
	if readLLength < lenToCompare:
		extraStart = lenToCompare - readLLength
		lenToCompare = readLLength
		# print extraStart
		# print_fullConsensus([], group.consensus)
		# print " "*offset, lread_L

	for i in xrange(lenToCompare):
		if len(group.consensus[i+len(seq_read_R)]) > 1 or \
				group.consensus[i+len(seq_read_R)].keys()[0] != \
				lread_L[i+lenCompared]:
			mismatches.add(i+len(seq_read_R))

	for m in group.mismatches:
		mismatches.add(m)
	for m in mis:
		mismatches.add(m)

	# if offset < 50:
	#	  print " "*offset, (seqs[read_L]+seqs[read_L+1])
	#	  print " "*-offset, seqs[read_R-1]+seqs[read_R]
	#	  print offset

	if len(mismatches) <= m2:
		group.mismatches = mismatches
	else:
		return False

	# Update consensus
	lenToUpdate = min(len(group.consensus)-offset, len(lread_L))
	for i in xrange(lenToUpdate):
		group.consensus[i+offset][lread_L[i]] = \
			  group.consensus[i+offset].get(lread_L[i], 0) + 1

	# Extend consensus to the right, if anything extra
	extraPart = len(lread_L) - len(group.consensus) + offset
	for i in xrange(extraPart):
		group.consensus.append({lread_L[-(extraPart-i)]:1})

	# if offset < 50:
	#	  print " "*offset, (seqs[read_L]+seqs[read_L+1])
	#	  print " "*-offset, seqs[read_R-1]+seqs[read_R]
	#	  print offset

	if group.leftReadsOffset > offset:
		group.leftReadsOffset = offset
		# print_fullConsensus([], group.consensus)
		# print " "*offset, lread_L
		# print offset
	elif group.maxLeftReadsOffset < offset:
		group.maxLeftReadsOffset = offset
		# print newLeftReadsOffset
		# print offset
		# print

	#group.mismatches = mismatches
	return True


def alignRightParts_secondPass(group, seqs):
	for pre_read_R in group.preRightParts:
		for offset in group.preRightParts[pre_read_R]:
			seq_pre_read_R = seqs[pre_read_R-1]+seqs[pre_read_R]
			shallowGroup = getShallowGroup(group, pre_read_R,
								seq_pre_read_R, offset, 0, 0, False)
			if shallowGroup:
				if pre_read_R in shallowGroup.preRightParts:
					shallowGroup.preRightParts[pre_read_R].add(offset)
				else:
					shallowGroup.preRightParts[pre_read_R] = set([offset])


def createNewGroup(group, next_read_R, seq_next_read_R, offset, mismatches):
	global c5
	c5 += 1

	newGroup = RightPartGroup(group.consensus)

	if offset > 0:
		j = offset
		k = 0
	else:
		j = 0
		k = -offset

	# if next_read_R == 1805:
	#	  print next_read_R
	#	  print "LOL-NEW"
	#	  print group.mismatches
	#	  print newGroup.mismatches
	#	  print mismatches
	#	  print

	# Update main-consensus
	for i in xrange(len(seq_next_read_R)-k):
		bp = seq_next_read_R[i+k]
		newGroup.consensus[i+j][bp] = newGroup.consensus[i+j].get(bp, 0) + 1

	# Update pre-consensus if any
	for i in xrange(k):
		bp = seq_next_read_R[i]
		newGroup.preConsensus.append({bp:1})

	global c2
	c2 += 1
	if next_read_R < secondSample:
		#newGroup.rightPartsN[next_read_R] = set([offset])
		newGroup.rightPartsN[next_read_R] = [offset]
	else:
		#newGroup.rightPartsD[next_read_R] = set([offset])
		newGroup.rightPartsD[next_read_R] = [offset]
	newGroup.mismatches = mismatches
	if len(mismatches) > M2:
		print mismatches



	group.rightPartGroups.append(newGroup)


def addToGroup(group, rightPartGroup, seqs, read_R, next_read_R, offset, mismatches, m2):
	global c4
	c4 += 1
	seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
	seq_read_R = seqs[read_R-1]+seqs[read_R]

	#all_mismatches = set([mis for mis in mismatches])
	all_mismatches = copy.deepcopy(mismatches)

	# Computes the length of the pre-consensus extension, if any
	lenPreConsensus = len(rightPartGroup.preConsensus)
	toExtend = -(lenPreConsensus + (group.readROffset -
				len(seqs[next_read_R-1]) + offset))

	if toExtend > 0:
		"""
		Case 1 - extending pre-consensus:
		 GAG TTATCATTGTGACTGGACAAAGTACG
		GGAG TTATCATTGTGACTGGACAAA
		"""
		beginning = lenPreConsensus
		j = 0
		l = 0
	elif offset > 0:
		"""
		Case 2 - no pre consensus access:
		GAG TTATCATTGTGACTGGACAAAGTACG
			   TCATTGTGACTGGACAAA
		"""
		beginning = 0
		j = offset
		toExtend = 0
	else:  # offset <= 0
		"""
		Case 3 - touching pre-consensus:
		GAG TTATCATTGTGACTGGACAAAGTACG
		 AG TTATCATTGTGACTGGACAAA
		"""
		# preconsensus already checked, start offset in read_k
		beginning = -offset
		# if read_k aligns in middle of consensus
		j = 0
		# Offset in preConsensus
		l = lenPreConsensus+offset
		# No extension
		toExtend = 0

	# Check consensus
	# start = beginning + toExtend
	# for i in xrange(group.leftReadsOffset-j):
	#	  if seq_read_R[i+j] != seq_next_read_R[i+start]:
	#		  all_mismatches.add(i+j)
	#		  if len(all_mismatches) > m2:
	#			  return True

	# Check consensus
	# start = beginning + toExtend
	# for i in xrange(group.leftReadsOffset-j):
	#	  if len(rightPartGroup.consensus[i+j]) > 1 or \
	#			 rightPartGroup.consensus[i+j].iterkeys().next() \
	#			 != seq_next_read_R[i+start]:
	#		  all_mismatches.add(i+j)
	#		  if len(all_mismatches) > m2:
	#			  return False


	# Check pre-consensus
	for i in xrange(beginning):
		if len(rightPartGroup.preConsensus[i+l]) > 1 or \
			   rightPartGroup.preConsensus[i+l].iterkeys().next() \
			   != seq_next_read_R[i+toExtend]:
			all_mismatches.add((i+l)-len(rightPartGroup.preConsensus))
			global c3
			c3 += 1
			if len(all_mismatches) > m2:
				# if next_read_R == 11035:
				#	  print lenPreConsensus
				#	  print "hej"
				#	  print offset
				#	  print toExtend
				#	  print beginning
				#	  print all_mismatches
				#	  print_fullConsensus(rightPartGroup.preConsensus, rightPartGroup.consensus)
				#	  print next_read_R
				# createNewGroup(group, next_read_R, seq_next_read_R, offset,
				#				 mismatches)
				return False


	# Add existing mismatches
	# for mis in group.mismatches:
	#	  all_mismatches.add(mis)
	for mis in rightPartGroup.mismatches:
		all_mismatches.add(mis)
		if len(all_mismatches) > m2:
			# createNewGroup(group, next_read_R, seq_next_read_R, offset,
			#				 mismatches)
			# global c3
			# c3 += 1
			# print read_R
			# print next_read_R
			# print " "*-offset,
			# print_fullConsensus(rightPartGroup.preConsensus, rightPartGroup.consensus)
			# print " "*offset, seq_next_read_R
			# print offset
			# print all_mismatches
			# print rightPartGroup.mismatches
			return False

	# Update pre-consensus
	for i in xrange(beginning):
		bp = seq_next_read_R[i+toExtend]
		rightPartGroup.preConsensus[i+l][bp] = \
			rightPartGroup.preConsensus[i+l].get(bp, 0)+1

	# Update main-consensus
	start = beginning + toExtend
	for i in xrange(len(seq_next_read_R)-start):
		bp = seq_next_read_R[i+start]
		rightPartGroup.consensus[i+j][bp] = \
			rightPartGroup.consensus[i+j].get(bp, 0) + 1

	# Extend pre-consensus if required
	if toExtend > 0:
		prePart = [{seq_next_read_R[i]:1} for i in xrange(toExtend)]
		rightPartGroup.preConsensus = prePart + rightPartGroup.preConsensus

	# Add read to group
	global c2
	c2 += 1
	# if next_read_R == 3449:
	#	  print next_read_R
	#	  print "LOL"
	#	  print group.mismatches
	#	  print rightPartGroup.mismatches
	#	  print
	if next_read_R < secondSample:
		if next_read_R in rightPartGroup.rightPartsN:
			#rightPartGroup.rightPartsN[next_read_R].add(offset)
			rightPartGroup.rightPartsN[next_read_R].append(offset)
		else:
			#rightPartGroup.rightPartsN[next_read_R] = set([offset])
			rightPartGroup.rightPartsN[next_read_R] = [offset]
	else:
		if next_read_R in rightPartGroup.rightPartsD:
			#rightPartGroup.rightPartsD[next_read_R].add(offset)
			rightPartGroup.rightPartsD[next_read_R].append(offset)
		else:
			#rightPartGroup.rightPartsD[next_read_R] = set([offset])
			rightPartGroup.rightPartsD[next_read_R] = [offset]

	rightPartGroup.mismatches = all_mismatches
	return True


def testRead(group, seqs, read_R, next_read_R, offset, m2, alignments, log):
	seq_read_R = seqs[read_R-1]+seqs[read_R]
	seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]

	# Check overlapping part
	leftROffset = group.leftReadsOffset
	lenToCompare = len(seq_next_read_R) - (leftROffset - offset)
	mismatches = set([mis for mis in group.mismatches])
	m1 = 0
	m = M1  # set to 0 or M1
	for i in xrange(lenToCompare):
		# if next_read_R == 1805:
		#	  print "hejhej"
		#	  print group.consensus[group.leftReadsOffset+i].keys()[0], seq_next_read_R[i+group.leftReadsOffset-offset]
		#	  print group.leftReadsOffset+i
		#	  print group.leftReadsOffset
		# print lenToCompare
		# print group.leftReadsOffset+i
		# print group.leftReadsOffset-offset+i
		# print group.consensus[i+group.leftReadsOffset].iterkeys().next()
		# print seq_next_read_R[i+group.leftReadsOffset-offset]

		if next_read_R < secondSample and (leftROffset+i) < len(seq_read_R):
			if seq_read_R[i+leftROffset] != \
				   seq_next_read_R[i+leftROffset-offset]:
				mismatches.add(leftROffset+i)
				m1 += 1
				if m1 > m:
					return alignments
		if len(group.consensus[leftROffset+i]) > 1 or \
				group.consensus[leftROffset+i].iterkeys().next() \
				!= seq_next_read_R[i+leftROffset-offset]:
			mismatches.add(leftROffset+i)
			global c3
			c3 += 1
			m1 += 1
			# if next_read_R == 1805:
			#	  #print " "*-offset,
			#	  print_fullConsensus([], group.consensus)
			#	  print " "*offset, seq_next_read_R
			#	  print offset
			#	  print "lololol"
			#	  print group.consensus[group.leftReadsOffset+i].keys()[0], seq_next_read_R[i+group.leftReadsOffset-offset]
			if leftROffset+i < len(seq_read_R):
				if m1 > M1:
					return alignments
			else:
				if m1 > M2:
					return alignments

	if len(mismatches) > m2:
		return alignments

	# print c3
	# sys.exit()
	# if next_read_R == 2097:
	#	  print " "*-offset, seq_read_R
	#	  print " "*offset, seq_next_read_R
	#	  print read_R
	#	  print next_read_R
	#	  print mismatches
	#	  print
	#	  # print "Passed"
	#	  print_fullConsensus([], group.consensus)
	#	  print group.maxLeftReadsOffset
	#	  print group.leftReadsOffset
	#	  print offset

	# Check rest against anchor
	if offset > 0:
		j = offset
		k = 0
	else:
		j = 0
		k = -offset
	for i in xrange(leftROffset-j):
		# if next_read_R == 2097:
		#	  print seq_read_R[i+j], seq_next_read_R[i+k]
		if seq_read_R[i+j] != seq_next_read_R[i+k]:
			global c3
			c3 += 1
			mismatches.add(i+j)
			if len(mismatches) > m2:
				return alignments

	added = False
	# Check if it fits into existing groups
	for rightPartGroup in group.rightPartGroups:
		if addToGroup(group, rightPartGroup, seqs, read_R, next_read_R,
					  offset, mismatches, m2):
			added = True
			break

	if not added:
		createNewGroup(group, next_read_R, seq_next_read_R, offset,
					   mismatches)

	return alignments + 1


def alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log):
	startOnePosOverlap = True
	for group in alignedGroups:
		group.checkedRightParts.add(read_R)
		# TODO: Should probably be removed as valid groups might be removed
		if len(group.leftPartsN)+len(group.leftPartsD) < 3:
			global c8
			c8 += 1
			continue
		for read_L in (group.leftPartsN.keys()+group.leftPartsD.keys()):
			if len(candidatePairs[read_L]) > maxCandMates:
				continue
			parts = candidatePairs[read_L]
			for next_read_R in parts:
				if next_read_R not in group.checkedRightParts:
					if startOnePosOverlap:
						# m = M2
						# if next_read_R < secondSample:
						#	  if M2 > 0 and len(group.mismatches) < M2:
						#		  m = M2-1
						alignments = 0
						seq_next_read_R = seqs[next_read_R-1] + \
										  seqs[next_read_R]
						#st = group.maxLeftReadsOffset
						# for offset in xrange(group.readROffset,
						#		  st-len(seq_next_read_R), -1):
						for offset in xrange(overlap, -overlap-1, -1):
						# offset = 0
						# for change in xrange(overlap*2+1):
						#	  if change % 2 == 0:
						#		  offset -= change
						#	  else:
						#		  offset += change
							# print "st:", st
							# print "r_offset:", group.readROffset
							# print_leftGroup(group, read_R, seqs, log)
							# print next_read_R
							# print seq_next_read_R
							# print offset
							# print st-len(seq_next_read_R)
							global numreadR
							numreadR += 1
							alignments = testRead(group, seqs, read_R,
												  next_read_R, offset, M2,
												  alignments, log)
							group.checkedRightParts.add(next_read_R)
							if alignments == maxAlignments:
								break
					else:
						for alignInfo in findAlignment(next_read_R, read_L,
									seqs, group.readROffset, 0, log):
							for offset2 in (group.leftPartsN[read_L] +
											group.leftPartsD[read_L]):
								# global numreadR
								# numreadR += 1
								offset1, mis, lenCompared = alignInfo
								offset = offset2 - offset1
								offset += group.readROffset - \
										  len(seqs[next_read_R-1])

								# mis = dismissRead(group, seqs, next_read_R,offset)
								#
								# if len(mis) > M2:
								#	  group.checkedRightParts.add(next_read_R)
								#	  continue
								#
								# added = False
								# for rightPartGroup in group.rightPartGroups:
								#	  added = fitsInGroup4(group, rightPartGroup,
								#		  seqs, read_R, next_read_R, offset,
								#		  mis, M2)
								#
								# # if not added:
								# if len(group.rightPartGroups) == 0:
								#	  seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
								#	  createNewGroup(group, next_read_R,
								#			  seq_next_read_R, offset, mis)

								testRead(group, seqs, read_R, next_read_R,
										 offset, M2, alignments, log)
								group.checkedRightParts.add(next_read_R)


		# Second pass to add possible missing pre-right-parts to other groups
		# alignRightParts_secondPass(group, seqs)
		# for shallowGroup in group.shallowGroups:
		#	  alignRightParts_secondPass(shallowGroup, seqs)

# GTCAA AGTTCAG
#		   TCAGAA TGCCC
#		  TTCAGAA TGCC
# 7 6 3 = 4
def newFindMutation(read_R, seqs, alignedGroups, log):
	numUsefulGroups = 0
	ancLen = len(seqs[read_R-1])+len(seqs[read_R])
	misInOverlap = 0
	for group in alignedGroups:
		for rightPartGroup in group.rightPartGroups:
			usefulMis = False
			for mis in rightPartGroup.mismatches:
				validPos = 0
				if mis < 0:
					for overlaps in rightPartGroup.preConsensus[mis].values():
						if overlaps >= requiredOverlaps:
							validPos += 1
				else:
					# if len(rightPartGroup.consensus[pos]) > 2:
					#	  print_alignedGroup(group, rightPartGroup, read_R,
					#						  seqs, log)
					for overlaps in rightPartGroup.consensus[mis].values():
						if overlaps >= requiredOverlaps:
							validPos += 1
				if validPos > 1:
					usefulMis = True
					if mis > group.leftReadsOffset and mis < ancLen:
						misInOverlap += 1
			if usefulMis and misInOverlap > 0:
				print_alignedGroup(group, rightPartGroup, read_R,
									seqs, log)
				numUsefulGroups += 1
				continue
	return numUsefulGroups


def oldFindMutation(read_R, seqs, alignedGroups, log):
	numUsefulGroups = 0
	anchorLen = len(seqs[read_R-1])+len(seqs[read_R])
	if log == None:
		log = open("usefulGroups/"+str(p_id)+".txt","w")
	for group in alignedGroups:
		for rightPartGroup in group.rightPartGroups:
			# Check if groups contains enough information to analyze
			numNormal = len(group.leftPartsN)+len(rightPartGroup.rightPartsN)
			if numNormal+1 < requiredOverlaps:
				continue
			numDiseas = len(group.leftPartsD)+len(rightPartGroup.rightPartsD)
			if numDiseas < requiredOverlaps:
				continue
			isUsefulGroup = False
			for mis in rightPartGroup.mismatches:
				if isUsefulGroup:
					break
				if mis < group.leftReadsOffset or mis > anchorLen-1:
					continue
				muts1 = findMutation(read_R, seqs, group.leftPartsN,
						 rightPartGroup.rightPartsN, mis, True, log)
				muts2 = findMutation(read_R, seqs, group.leftPartsD,
						 rightPartGroup.rightPartsD, mis, False, log)
				for mut1 in muts1:
					for mut2 in muts2:
						if mut1 != mut2 and not isUsefulGroup:
							read = seqs[read_R]
							print_compare(read, mis, mut1,
										  rightPartGroup.consensus, log)
							isUsefulGroup = True
							numUsefulGroups += 1
							print_alignedGroup(group, rightPartGroup, read_R,
													seqs, log)

	return numUsefulGroups


def findMutation(read_R, seqs, leftparts, rightparts, mutationsPos, first, log):
	# Test that the size of the group is big enough to infer any info
	# if (len(leftparts)+len(rightparts)) < 3:
	#	  return "Fail"
	count = 0
	validBPs = []
	#for mutationsPos in mismatches:
	# print "mutPos", mutationsPos
	#firstSampleBP = set()
	mutationBPsLeft = dict()
	mutationBPsRight = dict()
	# Get bp in anchor point
	if first and read_R < secondSample or \
			not first and read_R > secondSample:
		read = seqs[read_R-1]+seqs[read_R]
		if mutationsPos > 0 and mutationsPos < len(read):
			#firstSampleBP.add(read[mutationsPos])
			bp = read[mutationsPos]
			mutationBPsRight[bp] = mutationBPsRight.get(bp, 0) + 1
	# Get bp in leftparts
	for read_L in leftparts:
		for offset in leftparts[read_L]:
			offset += len(seqs[read_L])
			if mutationsPos < offset:
				continue
			read = seqs[read_L]+seqs[read_L+1]
			if mutationsPos > len(read)+offset-1:
				continue
			# print read[mutationsPos-offset-1]
			# print read[mutationsPos-offset]
			# print read[mutationsPos-offset+1]
			count += 1
			#firstSampleBP.add(read[mutationsPos-offset])
			bp = read[mutationsPos-offset]
			mutationBPsLeft[bp] = mutationBPsLeft.get(bp, 0) + 1
	# Get bp in rightparts
	for next_read_R in rightparts:
		for offset in rightparts[next_read_R]:
			if mutationsPos < offset:
				continue
			read = seqs[next_read_R-1]+seqs[next_read_R]
			if mutationsPos > len(read)+offset-1:
				continue
			# print mutationsPos-offset
			# print read[mutationsPos-offset-1]
			# print read[mutationsPos-offset]
			# print read[mutationsPos-offset+1]
			# print mutationsPos
			# print offset
			# print read_R
			# print read
			# print len(read)
			count += 1
			#firstSampleBP.add(read[mutationsPos-offset])
			bp = read[mutationsPos-offset]
			mutationBPsRight[bp] = mutationBPsRight.get(bp, 0) + 1
	for bp in mutationBPsRight:
		overlapsLeft = mutationBPsLeft.get(bp, 0)
		if overlapsLeft > 0:
			if mutationBPsRight[bp]+overlapsLeft >= requiredOverlaps:
				validBPs.append(bp)
	return validBPs


# ************************************************************************** #
#																			 #
#						   Sequence Align All Reads							 #
#																			 #
# ************************************************************************** #
def seqAlignAllReads(fasta_file, log):
	seqs = getAllReads(fasta_file, log)

	numParts = len(seqs)/2
	prog = 0
	tim = time.clock()
	for read_R in xrange(1,len(seqs),2):
		alignedGroups = []

		# Align left parts
		alignLeft(read_R, seqs, alignedGroups, log)

		# Align right parts
		alignRight(read_R, seqs, alignedGroups, log)

		# print_alignedGroups(alignedGroups, read_R, seqs, log)
		# sys.exit()

		prog += 1
		if prog % 500 == 0:
			logprint(log, False, "Processed", prog, "of", numParts,
					 "right parts in", (time.clock()-tim) / 60, "minutes")
			logprint(log, True, "Memory usage (in mb):",
					 memory_usage_resource())
			global c1
			logprint(log, False, "left parts aligned:", c1)
			c1 = 0
			global c2
			logprint(log, True, "right parts aligned:", c2)
			c2 = 0

	logprint(log, False, "Finished sequence alignment",
			 (time.clock() - tim) / 60, "minutes")
	logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
	logprint(log, False, "c1:", c1)
	logprint(log, False, "c2:", c2)
	logprint(log, False, "c3:", c3)
	logprint(log, False, "c4:", c4)
	logprint(log, False, "c5:", c5)
	logprint(log, False, "c6:", c6)
	logprint(log, False, "c7:", c7)
	logprint(log, False, "numReadL:", numreadL)
	logprint(log, False, "numReadR:", numreadR)


def alignLeft(read_R, seqs, alignedGroups, log):
	readROffset = len(seqs[read_R-1])
	for read_L in xrange(0,len(seqs),2):
		#print seqs[read_R]
		#print seqs[read_L]
		for alignInfo in findAlignment(read_R, read_L, seqs,
										 readROffset, M1, log):
			#offset += len(seqs[read_R-1])
			offset, mis, lenCompared = alignInfo
			global numreadL
			numreadL += 1
			#print offset
			newGroup = True
			for group in alignedGroups:
				if fitsInGroup(group, seqs, read_R, read_L, alignInfo, offset, M2):
					# Add read_L to group
					global c1
					c1 += 1
					if read_L in group.leftParts:
						group.leftParts[read_L].append(offset)
					else:
						group.leftParts[read_L] = [offset]

					# Fits in group, so don't create new group
					newGroup = False
					break
			# If read_L doesn't fit in any group, then create a new
			if newGroup:
				global c7
				c7 += 1
				group = AlignedGroup(read_R, readROffset, offset+readROffset)
				# Start of extra part - offsetExtraPart
				start = len(seqs[read_R]) - offset
				#group.consensusRight = seqs[read_L][start:] + seqs[read_L+1]

				# Add read_L to new the group
				group.leftParts[read_L] = [offset]
				group.leftReadsOffset = offset+group.readROffset
				group.mismatches = mis

				# Add anchor point to consensus
				for bp in ''.join((seqs[read_R-1], seqs[read_R])):
					group.consensus.append({bp:1})

				# Add overlapping part of read_L to consensus
				for index in xrange(start):
					i = index + group.readROffset + offset
					bp = seqs[read_L][index]
					group.consensus[i][bp] = group.consensus[i].get(bp, 0) + 1

				# Add the rest of read_L to consensus
				for bp in ''.join((seqs[read_L][start:], seqs[read_L+1])):
					group.consensus.append({bp:1})

				# Append new group to the other groups
				alignedGroups.append(group)


def alignRight(read_R, seqs, alignedGroups, log):
	for group in alignedGroups:
		for read_L in group.leftParts:
			for next_read_R in xrange(1,len(seqs),2):
				if read_R != next_read_R and \
							 next_read_R not in group.checkedRightParts:
					seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
					st = group.maxLeftReadsOffset
					for offset in xrange(st-len(seq_next_read_R)+10, st):
						global numreadR
						numreadR += 1

						testRead(group, seqs, read_R, next_read_R,
								 offset, M2, log)
						group.checkedRightParts.add(next_read_R)


# ************************************************************************** #
#																			 #
#									  Main									 #
#																			 #
# ************************************************************************** #
def main():
	"""
	Main method of the program
	"""
	totim = time.clock()

	# s1 = "TACTCGCAATATCAAGCCTTATGGTTTATTCGTTCCTATGATCAGGATGA"
	# s2 = "TCAAGCCTTATGGTTTATTCGTTCCTATGATCAGGATGATCATCAACTTC"
	# print globalAlignment(s1, s2, 38, True)
	# sys.exit()

	# Parse command line options
	fasta_file, normal_file, diseased_file, k, threshold, bands, rows, \
		 minhash_alg, seed, log_file, input_file, output_file, test \
		 = optionParse()

	# For testing only
	# input_file = "candidate_pairs_k_16_b_2_r_5_m_6.txt"

	# n (number of hash functions = length of minhash signature) is given by
	# n = bands * rows

	with open(log_file, 'w') as log:

		if input_file:
			candidatePairs = importCandidatePairs(input_file, log)
		else:
			# candidatePairs = runLSH(fasta_file, bands, rows, n, k, seed,
			#						  minhash_alg, log)
			multiProcessing = False
			if multiProcessing:
				p_size = bands
				# pool = Pool(processes=p_size)
				pool = None
				candidatePairs = runLSH(normal_file, diseased_file, bands,
					 				rows, k, seed, minhash_alg, test, log,
									multiProcessing, pool)

			else:
				candidatePairs = runLSH(normal_file, diseased_file, bands,
									rows, k, seed, minhash_alg, test, log,
									multiProcessing, None)

		if output_file:
			output = "candidate_pairs_k_"+str(k)+"_b_"+str(bands)+"_r_"+ \
					 str(rows)+"_m_"+str(minhash_alg)
			exportCandidatePairs(candidatePairs, output_file, log)

		if test == 0:
			"""
			Stop after the LSH step.
			"""
			print "Total time used:", time.clock() - totim / 60, "minutes"
			sys.exit()
		elif test == 1:
			"""
			Test S-curve distribution
			"""
			makeSPlot(normal_file, candidatePairs, k, bands, rows,
					  minhash_alg, log)
		elif test == 2:
			"""
			Test which pairs are found by LSH
			"""
			pairsFoundByLSH(normal_file, diseased_file, candidatePairs, k,
							bands, rows, log)
		elif test == 3:
			"""
			Test different k values
			"""
			if fasta_file:
				os.system("echo '"+str(k)+"\t"+str(candidatePairs[0])+"\t"+
						  str(candidatePairs[1])+"\t"+
						  str(candidatePairs[2])+"' >> "+fasta_file)
		elif test == 4:
			"""
			Count number of candidate pairs satisfying the naive sim
			"""
			naiveSim = 0.97
			total = 0
			satisfying = 0
			seqs = getAllReads(normal_file, log) + \
				   getAllReads(diseased_file, log)
			for doc1 in candidatePairs:
				for doc2 in candidatePairs[doc1]:
					naive = globalAlignment(seqs[doc1],seqs[doc2], 32)
					if naive >= naiveSim:
						satisfying += 1
					else:
						logprint(log, False, seqs[doc1], "\n"+seqs[doc2],
							 	 "\n"+str(naive), "\n")
					total += 1
			logprint(log, False, "Satisfying naive sim:", satisfying)
			logprint(log, False, "Total pairs:", total)
			logprint(log, False, "Ratio:", float(satisfying)/total)
		elif test == 5:
			p_size = bands
			# pool = Pool(processes=p_size)
			pool = None
			initMultiSeqAlign(normal_file, diseased_file, candidatePairs,
				 			  pool, p_size, log)
		else:
			# If no test to run on LSH, continue with constructing groups of aligned reads
            #findMutations(candidatePairs, normal_file, diseased_file, log)
			sequenceAlignment(candidatePairs, normal_file, diseased_file, log)

		logprint(log, True, "Total time used:", (time.clock() - totim) / 60,
				 "minutes")

		return (time.clock() - totim)



		### For Testing ###
		reads = getAllReads(fasta_file, log)

		findSimPairs = False
		if findSimPairs:
			findSimilarPairs(reads, candidatePairs, k, bands, rorws,
							 minhash_alg, log)
		else:
			validateJaccardSim(reads, candidatePairs, k, bands, rows,
							   minhash_alg, log)

		print "Total time used:", time.clock() - totim / 60, "minutes"


if __name__ == '__main__':
	"""
	Main method call
	"""
	main()

	sys.exit(0)
