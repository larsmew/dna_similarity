#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/12/2014"
__version__ = "$Revision: 2.0"

from optparse import OptionParser
from operator import itemgetter
from collections import Counter
from collections import deque
from itertools import izip
import sys
import time
import random
import json
import cPickle as pickle
import csv


c1 = 0
c2 = 0
M1 = 1
M2 = 1

# ************************************************************************** #
#                                                                            #
#                                   Helpers                                  #
#                                                                            #
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
                      help="set <VALUE> as similairy measure to use.")

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

    (options, args) = parser.parse_args()

    return options.fasta_file, options.k, options.threshold,\
        options.bands, options.rows, options.sim, options.m, \
        options.s, options.log, options.input, options.output


def memory_usage_resource():
    """
    Computes the ressource usage at a given time during runtime.
    Computes total amount used so far, so not the amount currently in use.
    """
    import resource
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
    #for element in output:
    #    log_output += str(element) + " "
    log_output = " ".join(map(str, output))
    print log_output
    log_output += "\n"
    log_file.write(log_output)
    if flush:
        log_file.flush()


def exportCandidatePairs(candidatePairs, output_file, log):
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
                f.write(str(sortedElements[-1])+"\n")

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
                elements = line.split("\t", 1)
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
    #     print id1
    #     print candidatePairs[id1]
    #     sys.exit()

    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
    return candidatePairs


# ************************************************************************** #
#                                                                            #
#                              Pre-computations                              #
#                                                                            #
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


def computeShinglesTable(fasta_file, k, log):
    """
    Computes a table for fast look-up of k-shingles and their corresponding
    position in the reads - when it was first encountered in the fasta file.
    """
    logprint(log, True, "Computing table of shingles positions...")
    tim = time.clock()
    with open(fasta_file, "rU") as fasta_file:
        shinglesPos = dict()
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


def computeShinglesSet(fasta_file, k, log):
    """
    Computes the set of all k-shingles (k-mers) in all reads.
    """
    logprint(log, True, "Computing set of all shingles...")
    tim = time.clock()
    with open(fasta_file, "rU") as fasta_file:
        shingles = set()
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
        return list(shingles)


def getAllReads(fasta_file, log):
    """
    Extract the reads (DNA sequences) from the given fasta file.
    """
    with open(fasta_file, "rU") as fasta_file:
        read = ""
        tim = time.clock()
        logprint(log, False, "Collecting reads...")
        reads = []
        seqs = 0
        for line in fasta_file:
            # If line starts with ">", which indicates end of a sequence, append it to list of reads
            if line.startswith(">"):
                if read != "":
                    seqs += 1
                    # Splits the string into two parts
                    leftpart = read[:len(read)/2]
                    rightpart = read[len(read)/2:]
                    reads.append(leftpart)
                    reads.append(rightpart)
                    read = ""
            # Concatenate multi-line sequences into one string
            else:
                read += line.strip().upper()
        if read != "":
            seqs += 1
            leftpart = read[:len(read)/2]
            rightpart = read[len(read)/2:]
            reads.append(leftpart)
            reads.append(rightpart)

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
    with open(fasta_file, "rU") as fasta_file:
        read = ""
        for line in fasta_file:
            # If line starts with ">", which indicates end of a sequence, append it to list of reads
            if line.startswith(">"):
                if read != "":
                    # Splits the string into two parts
                    leftpart = read[:len(read)/2]
                    yield leftpart
                    rightpart = read[len(read)/2:]
                    yield rightpart
                    read = ""
            # Concatenate multi-line sequences into one string
            else:
                read += line.strip().upper()
        if read != "":
            leftpart = read[:len(read)/2]
            yield leftpart
            rightpart = read[len(read)/2:]
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
#                                                                            #
#                         Locality Sensitive Hashing                         #
#                                                                            #
# ************************************************************************** #
def runLSH(fasta_file, bands, rows, n, k, seed, minhash_alg, log):
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

        Type 1:
            1. First hash - pre-computed hash functions
            2. First hash - Ongoing hash function
        Type 2:
            3. Whole matrix - pre-computed hash functions
            4. Whole matrix - Ongoing hash function
        Type 3:
            5. Doc shingles - pre-computed hash functions
            6. Doc shingles - Ongoing hash function
    """
    if fasta_file:
        tim = time.clock()
        random.seed(seed)
        candidatePairs = dict()

        # Computes table of all k-shingles and their position
        if minhash_alg == 3 or minhash_alg > 5:
            shingles = computeShinglesTable(fasta_file, k, log)
        # Computes set of all k-shingles
        else:  # minhash alg 1, 2, 4 or 5
            shingles = computeShinglesSet(fasta_file, k, log)
        # Use Locality-Sensitive Hashing computing for each bands the buckets
        # with similar documents (reads) obtained by minhashing each read.
        for b in xrange(bands):
            buckets = dict()
            minhashing(fasta_file, shingles, buckets, k, rows,
                       minhash_alg, b, bands, log)
            lshBand(buckets, b, candidatePairs, log)

        # Convert sets to lists for memory effciency
        # for id1 in candidatePairs:
        #     candidatePairs[id1] = list(candidatePairs[id1])

        logprint(log, False, "\nNumber of unique candidate pairs",
                 sum(len(candidatePairs[i]) for i in candidatePairs)/2)
        logprint(log, False, "Finished LSH in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())

        return candidatePairs

    else:
        logprint(log, True, "ERROR: NO FASTA FILE OR IMPORT FILE PROVIDED")
        sys.exit()


def getDocShingles(dna, k):
    """
    Computes all shingles of size k in a document (dna sequence)
    """
    shingles = {dna[i:i+k] for i in xrange(len(dna)-k+1)}
    return shingles


def minhashing(fasta, shingles, buckets, k, rows, minhash_alg, bn, bs, log):
    tim = time.clock()
    logprint(log, True, "Minhashing...")

    # minhashAlg = { '1' : minhashing_alg1,
    #                '2' : minhashing_alg2,
    #                '3' : minhashing_alg3,
    #                '4' : minhashing_alg4,
    #                '5' : minhashing_alg5,
    #                '6' : minhashing_alg6,
    # }

    idx = 0
    if minhash_alg < 4:
        hashfuncs = computeHashFunctions(rows, shingles, log)
    else:
        numShingles = len(shingles)
        p = getPrime(numShingles)
        a = [random.randrange(numShingles) for i in xrange(rows)]
        b = [random.randrange(numShingles) for i in xrange(rows)]
    for part in getPartsFromFile(fasta, log):
        if minhash_alg == 1:
            minhashing_alg1(part, idx, shingles, buckets, k, rows, hashfuncs)
        elif minhash_alg == 2:
            minhashing_alg2(part, idx, shingles, buckets, k, rows, hashfuncs)
        elif minhash_alg == 3:
            minhashing_alg3(part, idx, shingles, buckets, k, rows, hashfuncs)
        elif minhash_alg == 4:
            minhashing_alg4(part, idx, shingles, buckets, k, rows, p, a, b)
        elif minhash_alg == 5:
            minhashing_alg5(part, idx, shingles, buckets, k, rows, p, a, b)
        else:  # Default to minhash alg 6
            minhashing_alg6(part, idx, shingles, buckets, k, rows, p, a, b)

        idx += 1

        if idx % 500000 == 0:
            logprint(log, True, "Band", bn+1, "of", str(bs)+":",
                     "Processed", idx, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    logprint(log, False, "Finished minhashing in",
             (time.clock() - tim) / 60, "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def lshBand(buckets, b, candidatePairs, log):
    tim = time.clock()
    logprint(log, True, "Running LSH and finding similar pairs...")
    numPairsUnique = 0
    b += 1
    for bucket in buckets:
        for i in xrange(len(buckets[bucket])):
            id1 = buckets[bucket][i]
            for j in xrange(i+1, len(buckets[bucket])):
                id2 = buckets[bucket][j]
                if id1 % 2 == 0 and id2 % 2 == 1:
                    if id1 + 1 != id2:
                        if id1 in candidatePairs:
                            candidatePairs[id1].add(id2)
                        else:
                            candidatePairs[id1] = set([id2])
                        if id2 in candidatePairs:
                            candidatePairs[id2].add(id1)
                        else:
                            candidatePairs[id2] = set([id1])
                        numPairsUnique += 1
                if id1 % 2 == 1 and id2 % 2 == 0:
                    if id1 - 1 != id2:
                        if id1 in candidatePairs:
                            candidatePairs[id1].add(id2)
                        else:
                            candidatePairs[id1] = set([id2])
                        if id2 in candidatePairs:
                            candidatePairs[id2].add(id1)
                        else:
                            candidatePairs[id2] = set([id1])
                        numPairsUnique += 1

    logprint(log, True, "Number of buckets in band", str(b)+":", len(buckets))
    numPairs = 0
    for bucket in buckets:
        inBucket = buckets[bucket]
        numPairs += len(inBucket) * (len(inBucket)-1) / 2
    logprint(log, False, "Number of candidate pairs in band", str(b)+":",
             numPairs)
    logprint(log, True, "Number of unique candidate pairs in band",
             str(b)+":", numPairsUnique)

        # print "Finished LSH for band", b, "in", (time.clock() - tim) / 60, \
        #       "minutes"

    return None


# ************************************************************************** #
#                                                                            #
#                            Minhashing algorithms                           #
#                                                                            #
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


def minhashing_alg3(dna, idx, shinglePos, buckets, k, rows, hashfuncs):
    """
    Uses pre-computed hashFuncs and table to find original shingle position,
    then find new shingle with smallest position in hash function.
    """
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    signature = []
    docShingles = getDocShingles(dna, k)
    numShingles = len(shinglePos)
    for h in hashfuncs:
        minVal = numShingles+1
        for shingle in docShingles:
            pos = shinglePos[shingle]
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
    Runs through each hash (whole matrix) and saves smallest value which
    is in doc
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


def minhashing_alg6(dna, idx, shingles, buckets, k, rows, p, a, b):
    """
    DEFAULT MINHASH ALGORITHM
    Uses hash functions in the form ((a*pos+b) mod p) mod N,
    where a and b random integers and p is prime and p > N.
    Computes original position of shingle by finding all shingles and
    enumerates them, then store them in a table for fast look up.
    Tale is called shingles.
    """
    # Find signature for each document
    signature = []
    docShingles = getDocShingles(dna, k)
    numShingles = len(shingles)
    for i in xrange(rows):
        minVal = numShingles+1
        for shingle in docShingles:
            pos = shingles[shingle]
            val = ((a[i]*pos+b[i]) % p) % numShingles
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
#                                                                            #
#                             Similarity checkers                            #
#                                                                            #
# ************************************************************************** #
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


def globalAlignment(doc1, doc2, extraInfo=False):
    """
    Aligning sequences by using a sliding window approach.
    Returns the best score (matches / seqlength) between the two sequences.
    """
    start = 0
    if extraInfo:
        bestScore = (0, 0, 0)
    else:
        bestScore = 0
    seqLength = len(doc1.dna)-start
    while seqLength > 24:
        # print seqLength, bestScore[1]
        matches = 0
        for i in xrange(seqLength):
            # print len(doc1.dna)-start
            if doc1.dna[i] == doc2.dna[i+start]:
                matches += 1
        # print bestScore
        score = matches / float(seqLength)
        if extraInfo:
            if score > bestScore[0]:
                # print score, bestScore[0]
                # print seqLength, matches, bestScore[1]
                bestScore = (score, matches, seqLength)
                if bestScore[0] == 1.0:
                    return bestScore
        else:
            if score > bestScore:
                bestScore = score
                if bestScore == 1.0:
                    return bestScore
        start += 1
        seqLength = len(doc1.dna)-start
    return bestScore


def jaccardSim(doc1, doc2, k, sim_measure="standard"):
    """
    Computing the jaccard similarity.
    Option to use jaccard bag similarity or standard jaccard similarity.
    """
    # ## Bag jaccard sim ## #
    if sim_measure == "bag":
        counterA = doc1.counterShingles
        counterB = doc2.counterShingles
        intersection = sum((counterA & counterB).values())
        if intersection == 0:
            return 0
        union = len(doc1.dna) + len(doc2.dna) - intersection
        return float(intersection) / union

    # ## standard jaccard sim ## #
    if sim_measure == "standard":
        shingles1 = getDocShingles(doc1, k)
        shingles2 = getDocShingles(doc2, k)
        intersection = shingles1.intersection(shingles2)
        if len(intersection) == 0:
            return 0
        union = shingles1.union(shingles2)
        return float(len(intersection)) / len(union)


def validateJaccardSim(reads, candidatePairs, k, b, r, alg, log):
    """
    Method for used for obtaining information about the validity of the pairs
    found by LSH
    """
    ### Used for saving each individual pairs similarity
    # f1 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
    #           "_m"+str(alg)+"_all.txt", 'w')
    # f2 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
    #           "_m"+str(alg)+"_lsh.txt", 'w')
    # f3 = open("candidate_pairs_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
    #           "_m"+str(alg)+"_rest.txt", 'w')

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
    for doc1 in xrange(numReads):
        start = 1
        if doc1 % 2 == 0:
            start = 3
        for doc2 in xrange(doc1+start, numReads, 2):
            dna1 = reads[doc1]
            dna2 = reads[doc2]
            jaccard = jaccardSim(dna1, dna2, k)

            ### Saves the similarity for each pair, and for pairs found by LSH
            # if jaccard > 0:
            #     f1.write(str(pairNum) + " " + str(jaccard) + "\n")
            #     rest = True
            #     if doc1 in candidatePairs:
            #         if doc2 in candidatePairs[doc1]:
            #             f2.write(str(pairNum) + " " + str(jaccard) + "\n")
            #             rest = False
            #             count += 1
            #     if rest:
            #         f3.write(str(pairNum) + " " + str(jaccard) + "\n")
            #         count2 += 1
            #     pairNum += 1

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
                logprint(log, True, "Processed", process, "pairs in time:",
                      (time.clock() - timer), "Found", pairNum, "cand. pairs")

    f = open("s_shape_info_b"+str(b)+"_r"+str(r)+"_k"+str(k)+
             "_m"+str(alg)+"_small_581steps.txt", 'w')
    for i in xrange(len(allPairs)):
        if allPairs[i] == 0:
            f.write(str(0)+" "+str(0)+" "+str(0)+str(possibleSims[i])+"\n")
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
    #     print doc1.id, doc2.id
    #     print sig1
    #     print sig2
    #     print lol
    #     print math.sqrt(lol)
    # else:
    #     sys.exit(0)

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
    #     print doc1.id, doc2.id
    #     print sig1
    #     print sig2
    #     print seq
    #     print
    # else:
    #     sys.exit(0)
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
    #               xrange((len(seqA)+1))]
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
    #     print row

    # sys.exit(0)


# ************************************************************************** #
#                                                                            #
#                              Sequence Alignment                            #
#                                                                            #
# ************************************************************************** #
class AlignedGroup(object):
    """
    
    """
    #consensus = deque()
    consensus = []
    readROffset = 0
    consensusMain = 0  # read_R
    # consensusLeft = ""
    # consensusRight = ""
    leftParts = dict()  # reads in group coming from left part
    rightParts = dict()  # reads in group coming from right part
    mismatches = 0
    checked = False
    #alignMatrix = dict()  # dictionaty of offsets

    # Initializer
    def __init__(self, read_R, offset):
        #self.consensus = deque()
        self.consensus = []
        self.readROffset = offset
        self.consensusMain = read_R
        # self.consensusLeft = ""
        # self.consensusRight = ""
        self.leftParts = dict()
        self.rightParts = dict()
        self.mismatches = 0
        self.checked = False
        #self.alignMatrix = dict()


def sequenceAlignment(candidatePairs, fasta_file, log):
    seqs = getAllReads(fasta_file, log)
    #alignMatrix = dict()
    #alignedGroups = dict()

    doPrint = False
    numParts = len(candidatePairs) / 2
    proc = 0
    newOffset = 0
    tim = time.clock()
    for read_R in candidatePairs:
        if read_R % 2 == 1:
        #if read_R == 275:
            #consensus = []
            alignedGroups = dict()

            # Align left parts
            alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log)

            # print "", ''.join(consensus.keys()[0] for consensus in alignedGroups[read_R][0].consensus)
            # print "", ''.join(str(consensus.values()[0])+" " for consensus in alignedGroups[read_R][0].consensus), "\n"

            # Align right parts
            alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log)

            if doPrint:
                for group in alignedGroups[read_R]:
                    if group.mismatches > 0:
                        logprint(log, False, "\nConsensus:")
                        logprint(log, False, "", ''.join(consensus.keys()[0] 
                                 for consensus in group.consensus))
                        newOffset = len(group.consensus[:group.readROffset])
                        for read_L in group.leftParts:
                            for offset in group.leftParts[read_L]:
                                logprint(log, False, " " * (newOffset +
                                     offset), seqs[read_L]+""+seqs[read_L+1])
                        for read_R2 in group.rightParts:
                            for offset in group.rightParts[read_R2]:
                                logprint(log, False, " " * (offset + 
                                         group.readROffset - len(seqs[read_R-1])),
                                         seqs[read_R2-1]+""+seqs[read_R2])
                        logprint(log, False, "read_R:", read_R)
                        logprint(log, False, "Left parts:",
                                 sorted(list(group.leftParts)))
                        logprint(log, False, "Number of left parts:",
                                 len(group.leftParts))
                        logprint(log, False, "Right parts:",
                                 sorted(list(group.rightParts)))
                        logprint(log, False, "Number of right parts:",
                                 len(group.rightParts))
                        # logprint(log, False, len(group.rightParts), "\n")

            # logprint(log, True, "read_R:", read_R)
            # for group in alignedGroups[read_R]:
            #     totalSize = len(group.leftParts) + len(group.rightParts)
            #     if totalSize > 9:
            #         logprint(log, False, "Total reads:", totalSize)
            #         logprint(log, False, len(group.leftParts))
            #         logprint(log, False, len(group.rightParts), "\n")

            # sys.exit()

            proc += 1
            if proc % 500000 == 0:
                logprint(log, False, "Processed", proc, "of", numParts,
                         "right parts in", (time.clock()-tim) / 60, "minutes")
                logprint(log, True, "Memory usage (in mb):",
                         memory_usage_resource())

    logprint(log, False, "Finished sequence alignment",
             (time.clock() - tim) / 60, "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
    logprint(log, False, "c1:", c1)
    logprint(log, False, "c2:", c2)


def getSequences(fasta_file, ids, log):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    with open(fasta_file, "rU") as fasta_file:
        read = ""
        tim = time.clock()
        #logprint(log, False, "Collecting sequences...")
        reads = []
        seqs = 0
        id = 0
        idx = 0
        for line in fasta_file:
            # If line starts with ">", which indicates end of a sequence, append it to list of reads
            if line.startswith(">"):
                if read != "":
                    # Splits the string into two parts
                    if idx == ids[id]:
                        leftpart = read[:len(read)/2]
                        reads.append(leftpart)
                        seqs += 1
                        id += 1
                        if id == len(ids):
                            read = ""
                            break
                    idx += 1
                    if idx == ids[id]:
                        rightpart = read[len(read)/2:]
                        reads.append(rightpart)
                        seqs += 1
                        id += 1
                        if id == len(ids):
                            read = ""
                            break
                    idx += 1
                    read = ""
            # Concatenate multi-line sequences into one string
            else:
                read += line.strip().upper()
        if read != "":
            if idx == ids[id]:
                leftpart = read[:len(read)/2]
                reads.append(leftpart)
                seqs += 1
                id += 1
            idx += 1

            if id != len(ids) and idx == ids[id]:
                rightpart = read[len(read)/2:]
                reads.append(rightpart)
                seqs += 1

        # logprint(log, False, "Finished reading in",
        #          (time.clock() - tim) / 60, "minutes")
        # logprint(log, False, "Found", seqs, "sequences in fasta file")
        # logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
        return reads


def alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log):
    alignedGroups[read_R] = []
    for read_L in candidatePairs[read_R]:
        for offset in findAlignment(seqs[read_R], seqs[read_L], M1, log):
            newGroup = True
            for group in alignedGroups[read_R]:
                if fitsInGroup(group, seqs, read_R, read_L, offset, M2):
                    # Add read_L to group
                    if read_L in group.leftParts:
                        group.leftParts[read_L].append(offset)
                    else:
                        group.leftParts[read_L] = [offset]
                    
                    # Don't create new group
                    newGroup = False
                    break
            # If read_L doesn't fit in any group, then create a new
            if newGroup:
                group = AlignedGroup(read_R, len(seqs[read_R-1]))
                # Start of extra part - offsetExtraPart
                start = len(seqs[read_R]) - offset
                #group.consensusRight = seqs[read_L][start:] + seqs[read_L+1]
                
                # Add read_L to new the group
                group.leftParts[read_L] = [offset]
                
                # Add anchor point to consensus
                for bp in ''.join((seqs[read_R-1], seqs[read_R])):
                    group.consensus.append({bp:1})
                # Add overlapping part of read_L to consensus
                for index, bp in enumerate(seqs[read_L][:start]):
                    i = index+group.readROffset+offset
                    group.consensus[i][bp] = group.consensus[i].get(bp, 0) + 1
                # Add the rest of read_L to consensus
                for bp in ''.join((seqs[read_L][start:], seqs[read_L+1])):
                    group.consensus.append({bp:1})
                
                # Append new group to the other groups
                alignedGroups[read_R].append(group)


                # key, = group.consensus[0]
                # print key
                # consensus = []
                # for pos in group.consensus:
                #     for bp in pos:
                #         consensus.append(bp)
                # print ''.join(consensus)
                # print ''.join((seqs[read_R-1], seqs[read_R],
                #                seqs[read_L][start:], seqs[read_L+1]))
                # sys.exit()
                
                # lread_R = seqs[read_R]
                # lread_L = seqs[read_L]
                # print "", lread_R
                # print " "*offset, lread_L
                # print " "*(offset), lread_L[len(lread_R)-offset:]
                # #print " "*(offset), lread_L[-extraPart:]
                # consensusRight = alignedGroups[read_R][-1].consensusRight
                # print " "*(len(lread_R)), consensusRight
                # lread_L = lread_L+seqs[read_L+1]
                # print " "*(offset), lread_L
                # print " "*(len(lread_R)), lread_L[len(lread_R)-offset:]
                # sys.exit()


def findAlignment(read_R, read_L, m1, log):
    offset = 0
    doPrint = False
    if len(read_R) > len(read_L):
        offset = len(read_R) - len(read_L)
        lengthToCompare = len(read_L)
    else:
        lengthToCompare = len(read_R)
    while lengthToCompare > 10:
        mismatches = 0
        for i in xrange(lengthToCompare):
            if read_R[i+offset] != read_L[i]:
                mismatches += 1
                if mismatches > m1:
                    break
        if mismatches <= m1:
            if doPrint:
                print "offset:", offset
                print read_R
                if offset == 0:
                    print read_L
                else:
                    print " "*(offset-1), read_L
            yield offset
        offset += 1
        lengthToCompare -= 1
    #yield -1


def fitsInGroup(group, seqs, read_R, read_L, offset, m2):
    lread_R = seqs[read_R]
    lread_L = seqs[read_L]+seqs[read_L+1]
    extraPart = len(lread_L)-len(lread_R)+offset
    mismatches = 0
    # print "", ''.join((consensus.keys()[0] for consensus in group.consensus))
    # print " "*group.readROffset,lread_R
    # print " "*(offset+group.readROffset), lread_L
    # print " "*(len(lread_R)), lread_L[len(lread_R)-offset:]
    # print " "*(len(lread_R)), lread_L[-extraPart:]
    # print extraPart
    # print " "*(len(lread_R)), group.consensusRight
    # print len(group.consensusRight)
    # sys.exit()
    #for i in xrange(min(extraPart, len(group.consensusRight))):
    #    if lread_L[i+len(lread_R)-offset] != group.consensusRight[i]:
    offset = group.readROffset + offset
    for pos, bp in izip(group.consensus[offset:], lread_L):
        #print pos, bp
        if len(pos) > 2 or pos.keys()[0] != bp:
            mismatches += 1
            if mismatches > m2:
                return False
    
    # Must fit in group if here, so update consensus
    for pos, bp in izip(group.consensus[offset:], lread_L):
        pos[bp] = pos.get(bp, 0) + 1
    
    # Add extra part of read_L to consensus, if any
    extraPart = len(lread_L) - len(group.consensus[offset:])
    if extraPart > 0:
        for bp in lread_L[-extraPart:]:
            group.consensus.append({bp:1})
            # print "", ''.join(consensus.keys()[0] for consensus in group.consensus)
            # print "", ''.join(str(consensus.values()[0]) for consensus in group.consensus)
    group.mismatches = mismatches
    return True


def fitsInGroup2(group, seqs, read_R, next_read_R, offset, offset2, m2):
    ### Check extra part at the right side, however is checked by
    ### findAlignment already, so probably not needed
    lread_R = seqs[read_R]
    lread_R2 = seqs[next_read_R]
    extraPart = offset
    mismatches = 0
    for i in xrange(min(extraPart, len(group.consensusRight))):
        if lread_R2[i+len(lread_R)-offset] != group.consensusRight[i]:
            mismatches += 1
            if mismatches > m2:
                return False

    ### Check extra part at the left side, first for anchor point
    lread_R = seqs[read_R]
    lread_R_l = seqs[read_R-1]
    lread_R2 = seqs[next_read_R-1]+seqs[next_read_R]
    # print " "*(len(seqs[next_read_R-1])), lread_R
    # print "", seqs[read_R-1]+lread_R
    # print " "*(offset), lread_R2
    # print offset2
    for i in xrange(offset2):
        # print lread_R2[i+len(seqs[read_R-1])-offset]
        # print lread_R[i]
        if lread_R2[i+len(seqs[read_R-1])-offset] != lread_R[i]:
            mismatches += 1
            if mismatches > m2:
                return False

    ### Check extra part at the left side, then for left part of anchor point
    # print "", lread_R_l
    # print " "*offset, lread_R2
    for i in xrange(len(lread_R_l)-offset):
        # print lread_R2[i]
        # print lread_R_l[i+offset]
        if lread_R2[i] != lread_R_l[i+offset]:
            mismatches += 1
            if mismatches > m2:
                return False

    return True


def fitsInGroup3(group, seqs, read_R, next_read_R, offset, offset2, m2):
    doPrint = False
    lread_R = seqs[read_R]
    lread_R_l = seqs[read_R-1]
    lread_R2 = seqs[next_read_R-1]+seqs[next_read_R]
    offsetPos = abs(offset)

    # print " "*(len(lread_R_l)+offsetPos), lread_R
    # print " "*offsetPos, lread_R_l
    # print "", seqs[next_read_R-1]+" "+seqs[next_read_R]

    mismatches = 0
    for i in xrange(offset2):
        # print lread_R2[i+len(seqs[next_read_R-1])+offsetPos]
        # print lread_R[i]
        if lread_R2[i+len(seqs[next_read_R-1])+offsetPos] != lread_R[i]:
            mismatches += 1
            if mismatches > m2:
                return False

    for i in xrange(len(lread_R_l)):
        # print lread_R2[i+offsetPos]
        # print lread_R_l[i]
        if lread_R2[i+offsetPos] != lread_R_l[i]:
            mismatches += 1
            if mismatches > m2:
                return False

    # group.consensusLeft = "GGGGAAG"
    # lread_R2 = "GAA"+lread_R2
    # offsetPos += 3
    # print " ", group.consensusLeft+lread_R_l
    # print " "*(len(group.consensusLeft)-offsetPos), lread_R2
    if offsetPos > len(group.consensusLeft):
        preOffset = offsetPos - len(group.consensusLeft)
        for i in xrange(len(group.consensusLeft)):
            # print lread_R2[i+preOffset]
            # print group.consensusLeft[i]
            if lread_R2[i+preOffset] != group.consensusLeft[i]:
                mismatches += 1
                if mismatches > m2:
                    if doPrint:
                        print "ALIGNS NICELY UNTIL PRE-PART! - preOffset bigger"
                        print " "*preOffset, group.consensusLeft+seqs[read_R-1]+seqs[read_R]+group.consensusRight
                        print "", lread_R2
                    #sys.exit()
                    global c1
                    c1 += 1
                    return False

        group.consensusLeft = lread_R2[:offsetPos]

    else:
        preOffset = len(group.consensusLeft) - offsetPos
        for i in xrange(offsetPos):
            # print lread_R2[i]
            # print group.consensusLeft[i+preOffset]
            if lread_R2[i] != group.consensusLeft[i+preOffset]:
                mismatches += 1
                if mismatches > m2:
                    if doPrint:
                        print "ALIGNS NICELY UNTIL PRE-PART! - groupOffset bigger"
                        print "", group.consensusLeft+seqs[read_R-1]+seqs[read_R]+group.consensusRight
                        print " "*preOffset, lread_R2
                    #sys.exit()
                    global c2
                    c2 += 1
                    return False
    return True


def fitsInGroup4(group, seqs, read_R, next_read_R, offset, offset2, m2):
    seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
    mismatches = 0
    # print " "*abs(offset), ''.join(consensus.keys()[0] for consensus in group.consensus)
    # print " "*abs(offset), seqs[read_R-1]+seqs[read_R]
    # print "", seq_next_read_R
    toExtend = len(seq_next_read_R) - len(group.consensus[:group.readROffset+len(seqs[read_R])+offset])
    # print seq_next_read_R
    # slicedCons = group.consensus[:group.readROffset+len(seqs[read_R])+offset]
    # print "", ''.join(i.keys()[0] for i in slicedCons)
    # print len(seq_next_read_R)
    # print len(slicedCons)
    # print toExtend
    # print len(group.consensus)-(group.readROffset+len(seqs[read_R])+offset)
    # print len(group.consensus)-\
    #       (len(group.consensus)-(group.readROffset+len(seqs[read_R])+offset))
    # print len(group.consensus)
    # sys.exit()
    if toExtend > 0:
        # global c1
        # c1 += 1
        #print " "*toExtend, seq_next_read_R[toExtend:]
        
        for pos, bp in izip(group.consensus, seq_next_read_R[toExtend:]):
            if len(pos) > 2 or pos.keys()[0] != bp:
                mismatches += 1
                if mismatches > m2:
                    return False
        
        for pos, bp in izip(group.consensus, seq_next_read_R[toExtend:]):
            pos[bp] = pos.get(bp, 0) + 1
            
        # Extend consensus with pre-part
        group.readROffset += toExtend
        prePart = [{bp:1} for bp in seq_next_read_R[:toExtend]]
        group.consensus = prePart + group.consensus
        # print "", ''.join(consensus.keys()[0] for consensus in group.consensus)
        # print "", ''.join(str(consensus.values()[0])+" " for consensus in group.consensus), "\n"
    else:
        # global c2
        # c2 += 1
        # print "", ''.join(consensus.keys()[0] for consensus in group.consensus)
        # print "", ''.join(str(consensus.values()[0])+" " for consensus in group.consensus)
        # print " "*(group.readROffset-len(seqs[read_R-1])), seqs[read_R-1]+seqs[read_R]
        # print " "*(group.readROffset-len(seqs[read_R-1])+offset), seq_next_read_R
        offset += group.readROffset-len(seqs[read_R-1])
        for pos, bp in izip(group.consensus[offset:], seq_next_read_R):
            if len(pos) > 2 or pos.keys()[0] != bp:
                mismatches += 1
                if mismatches > m2:
                    return False

        for pos, bp in izip(group.consensus[offset:], seq_next_read_R):
            pos[bp] = pos.get(bp, 0) + 1
        group.mismatches += mismatches
    return True


def alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log):
    newOffset = 0
    for group in alignedGroups[read_R]:
        for read_L in group.leftParts:
            for next_read_R in candidatePairs[read_L]:
                if read_R != next_read_R and \
                             next_read_R not in group.rightParts:
                    for offset1 in findAlignment(seqs[next_read_R], seqs[read_L], 0, log):
                        for offset2 in group.leftParts[read_L]:
                            offset = offset2 - offset1
                            # if offset > 0:
                            #
                            #     if fitsInGroup2(group, seqs, read_R, next_read_R, offset, offset2, 0):
                            #         if next_read_R in group.rightParts:
                            #             group.rightParts[next_read_R].append(offset)
                            #         else:
                            #             group.rightParts[next_read_R] = [offset]
                            #
                            # elif offset < 0:
                            #
                            #     if fitsInGroup3(group, seqs, read_R, next_read_R, offset, offset2, 0):
                            #         if next_read_R in group.rightParts:
                            #             group.rightParts[next_read_R].append(offset)
                            #         else:
                            #             group.rightParts[next_read_R] = [offset]

                            if fitsInGroup4(group, seqs, read_R, next_read_R, offset, offset2, M2):
                                group.rightParts[next_read_R] = [offset]
    return newOffset

# GTCAA AGTTCAG
#          TCAGAA TGCCC
#         TTCAGAA TGCC
# 7 6 3 = 4

# ************************************************************************** #
#                                                                            #
#                                     Main                                   #
#                                                                            #
# ************************************************************************** #
def main():
    """
    Main method of the program
    """
    totim = time.clock()

    # Parse command line options
    fasta_file, k, threshold, bands, rows, sim_measure, minhash_alg, \
        seed, log_file, input_file, output_file = optionParse()

    # Temp paramters
    # input_file = "candidate_pairs_k_16_b_2_r_5_m_6"
    # input_file = "candidate_pairs_k_16_b_1_r_32_m_6"

    # n (number of hash functions = length of minhash signature) is given by
    n = bands * rows

    with open(log_file, 'w') as log:

        if input_file:
            candidatePairs = importCandidatePairs(input_file, log)
        else:
            candidatePairs = runLSH(fasta_file, bands, rows, n, k, seed,
                                    minhash_alg, log)
        if output_file:
            output = "candidate_pairs_k_"+str(k)+"_b_"+str(bands)+"_r_"+ \
                     str(rows)+"_m_"+str(minhash_alg)
            exportCandidatePairs(candidatePairs, output_file, log)

        sequenceAlignment(candidatePairs, fasta_file, log)

        logprint(log, True, "Total time used:", (time.clock() - totim) / 60,
                 "minutes")

        sys.exit()
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

    main()

    sys.exit(0)
