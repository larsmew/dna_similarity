#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/12/2014"
__version__ = "$Revision: 2.0"

from optparse import OptionParser
from operator import itemgetter
from collections import Counter
import sys
import time
import random
import copy
import json
import cPickle as pickle
import csv


c1 = 0
c2 = 0
c3 = 0
c4 = 0
c5 = 0
c6 = 0
c7 = 0
numreadL = 0
numreadR = 0
M1 = 1
M2 = 1

old = False

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
    #signature = array.array('l')
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
    consensus = []
    #preConsensus = []
    readROffset = 0
    consensusMain = 0  # read_R
    leftParts = dict()  # reads in group coming from left part
    leftReadsOffset = 0
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
        self.leftParts = dict()
        self.leftReadsOffset = leftReadsOffset
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
    rightParts = dict()
    mismatches = set()

    # Initializer
    def __init__(self, consensus):
        self.consensus = copy.deepcopy(consensus)
        self.preConsensus = []
        self.rightParts = dict()
        self.mismatches = set()


def print_fullConsensus(preconsensus, consensus, log=None):
    alphabetSize = 4
    for i in xrange(alphabetSize):
        consensusString = ""
        for j in xrange(len(preconsensus)):
            if i < len(preconsensus[j]):
                consensusString += preconsensus[j].keys()[i]
            else:
                consensusString += " "
        if consensusString != "":
            consensusString += " "
        for j in xrange(len(consensus)):
            if i < len(consensus[j]):
                consensusString += consensus[j].keys()[i]
            else:
                consensusString += " "
        if consensusString.strip() != "":
            print "", consensusString
            if log:
                log.write(" "+consensusString+"\n")
        else:
            break


def print_alignedGroups(groups, read_R, seqs, log):
        #for group in alignedGroups:
        for group in groups:
            if len(group.rightPartGroups) > 0:
                for rightPartGroup in group.rightPartGroups:
                    if len(rightPartGroup.mismatches) > -1:
                        logprint(log, False, "\nread_R:", read_R)
                        logprint(log, False, "Consensus:")
                        # logprint(log, False, "", ''.join(consensus.keys()[0]
                        #         for consensus in group.preConsensus) + #" " +
                        #         ''.join(consensus.keys()[0] for consensus
                        #         in group.consensus))
                        print_fullConsensus(rightPartGroup.preConsensus,
                                            rightPartGroup.consensus, log)
                        lenPre = len(rightPartGroup.preConsensus)
                        newOffset = lenPre + group.readROffset
                        for read_L in group.leftParts:
                            for offset in group.leftParts[read_L]:
                                logprint(log, False, " " * (newOffset +
                                     offset), seqs[read_L]+""+seqs[read_L+1])
                        for read_R2 in rightPartGroup.rightParts:
                            for offset in rightPartGroup.rightParts[read_R2]:
                                logprint(log, False, " " * (offset + lenPre),
                                         seqs[read_R2-1]+""+seqs[read_R2])
                        logprint(log, False, "Left parts:",
                                 sorted(list(group.leftParts)))
                        logprint(log, False, "Number of left parts:",
                                 len(group.leftParts))
                        logprint(log, False, "Right parts:",
                                 sorted(rightPartGroup.rightParts))
                        logprint(log, False, "Number of right parts:",
                                 len(rightPartGroup.rightParts))
                        logprint(log, False, "mismatches:",
                                 list(rightPartGroup.mismatches))
                        # logprint(log, False, "Number of mismatches:",
                        #          len(rightPartGroup.mismatches))
                        # if len(group.mismatches) > 0:
                        #     sys.exit()
            # else:
            elif len(group.mismatches) > 0:
                logprint(log, False, "\nread_R:", read_R)
                logprint(log, False, "Consensus:")
                # logprint(log, False, "", ''.join(consensus.keys()[0]
                #          for consensus in group.consensus))
                print_fullConsensus([],
                                    group.consensus, log)
                for read_L in group.leftParts:
                    for offset in group.leftParts[read_L]:
                        logprint(log, False, " " *(offset +
                                 group.readROffset),
                                 seqs[read_L]+""+seqs[read_L+1])
                logprint(log, False, "mismatches:", group.mismatches)
                logprint(log, False, "Left parts:",
                         sorted(list(group.leftParts)))
                logprint(log, False, "Number of left parts:",
                         len(group.leftParts))


def sequenceAlignment(candidatePairs, fasta_file, log):
    seqs = getAllReads(fasta_file, log)

    numParts = len(candidatePairs) / 2
    prog = 0
    tim = time.clock()
    for read_R in candidatePairs:
        if read_R % 2 == 1:
        # if read_R == 577:
            alignedGroups = []

            # Align left parts
            alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log)

            # Align right parts
            alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log)

            print_alignedGroups(alignedGroups, read_R, seqs, log)
            sys.exit()

            prog += 1
            if prog % 500000 == 0:
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


def alignLeftParts(read_R, seqs, alignedGroups, candidatePairs, log):
    readROffset = len(seqs[read_R-1])
    for read_L in candidatePairs[read_R]:
        #print seqs[read_R]
        #print seqs[read_L]
        for alignInfo in findAlignment(seqs[read_R], seqs[read_L],
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


def findAlignment(read_R, read_L, readROffset, m1, log):
    offset = 0
    doPrint = False
    if len(read_R) > len(read_L):
        offset = len(read_R) - len(read_L)
        lengthToCompare = len(read_L)
    else:
        lengthToCompare = len(read_R)
    while lengthToCompare > 10:
        mismatches = set()
        if log == 2119:
            print "", read_R
            print " "*offset, read_L
        for i in xrange(lengthToCompare):
            if read_R[i+offset] != read_L[i]:
                mismatches.add(i+offset+readROffset)
                if len(mismatches) > m1:
                    break
        if len(mismatches) <= m1:
            if doPrint:
                print "offset:", offset
                print read_R
                if offset == 0:
                    print read_L
                else:
                    print " "*(offset-1), read_L
                print mismatches
            yield offset, mismatches, lengthToCompare
        offset += 1
        lengthToCompare -= 1
    #yield -1


def fitsInGroup(group, seqs, read_R, read_L, alignInfo, offset, m2):
    global c6
    c6 += 1
    lread_R = seqs[read_R]
    lread_L = seqs[read_L]+seqs[read_L+1]
    mismatches = 0
    offset += group.readROffset
    offset2, mis, lenCompared = alignInfo

    if old:
        # Check consensus for mismatches before alignment of read_L
        # for i in xrange(offset-group.readROffset):
        #     if len(group.consensus[i+group.readROffset]) > 1:
        #         mismatches += 1
        #         if mismatches > m2:
        #             return False

        lenToCompare = min(len(group.consensus)-offset, len(lread_L))
        for i in xrange(lenToCompare):
            if len(group.consensus[i+offset]) > 1 or \
               group.consensus[i+offset].keys()[0] != lread_L[i]:
                mismatches += 1
                if mismatches > m2:
                    return False

        # Check consensus for mismatches after alignment of read_L
        # for i in xrange(len(group.consensus) - lenToCompare - offset):
        #     if len(group.consensus[i + lenToCompare + offset]) > 1:
        #         mismatches += 1
        #         if mismatches > m2:
        #             return False
    else:
        mismatches = set()
        seq_read_R = seqs[read_R-1]+seqs[read_R]
        # lenToCompare = min(len(group.consensus)-offset, len(lread_L))
        # for i in xrange(lenToCompare):
        #     if len(group.consensus[i+offset]) > 1 or \
        #        group.consensus[i+offset].keys()[0] != lread_L[i]:
        #         mismatches.add(i+offset)

        lenToCompare = len(group.consensus)-len(seq_read_R)
        readLLength = len(lread_L) - lenCompared
        # print lenCompared
        # print readLLength
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

    newLeftReadsOffset = offset+group.readROffset
    if group.leftReadsOffset > newLeftReadsOffset:
        group.leftReadsOffset = newLeftReadsOffset
    # if group.leftReadsOffset > offset:
    #     group.leftReadsOffset = offset
        # print_fullConsensus([], group.consensus)
        # print " "*offset, lread_L
        # print offset

    #group.mismatches = mismatches
    return True


def fitsInGroup4(groupOrg, group, seqs, read_R, next_read_R, offset, mis, m2):
    global c4
    c4 += 1
    # print_fullConsensus(group.preConsensus, group.consensus)
    seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
    seq_read_R = seqs[read_R-1]+seqs[read_R]

    # Computes the length of the pre-consensus extension, if any
    lenPreConsensus = len(group.preConsensus)
    toExtend = -(lenPreConsensus + (groupOrg.readROffset -
                len(seqs[next_read_R-1]) + offset))

    # Computes real offset of next_read_R in consensus - only changes if reads
    # have varying lengths
    # offset += group.readROffset - len(seqs[next_read_R-1])

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

    if old:

        mismatches = 0

        # Check if next_read_R matches main-consensus
        start = beginning + toExtend
        for i in xrange(len(seq_next_read_R)-start):
            if len(group.consensus[i+j]) > 1 or group.\
                   consensus[i+j].keys()[0] != seq_next_read_R[i+start]:
                mismatches += 1
                if mismatches > m2:
                    return False

        # Check positions in pre-consensus before next_read_R for mismatches
        for i in xrange(lenPreConsensus-beginning):
            if len(group.consensus[i]) > 1:
                mismatches += 1
                if mismatches > m2:
                    return False

        # Check positions in consensus before next_read_R for mismatches
        for i in xrange(offset):
            if len(group.consensus[i]) > 1:
                mismatches += 1
                if mismatches > m2:
                    return False

        # Check if next_read_R matches pre-consensus
        for i in xrange(beginning):
            if len(group.preConsensus[i+l]) > 1 or group.\
                   preConsensus[i+l].keys()[0] != seq_next_read_R[i+toExtend]:
                mismatches += 1
                if mismatches > m2:
                    global c3
                    c3 += 1
                    return False
    else:
        localMismatches = set([mismatch for mismatch in mis])
        mismatches = mis
        for mismatch in groupOrg.mismatches:
            mismatches.add(mismatch)
        for mismatch in group.mismatches:
            mismatches.add(mismatch)


        # print "hej"
        # print mismatches

        # if offset < 0:
        #     print " "*(-offset), print_fullConsensus([],group.consensus)
        # else:
        #     print_fullConsensus([],group.consensus)
        # print " "*offset, seq_next_read_R

        start = beginning + toExtend
        for i in xrange(groupOrg.leftReadsOffset):
            if len(group.consensus[i+j]) > 1 or group.\
                   consensus[i+j].keys()[0] != seq_next_read_R[i+start]:
                mismatches.add(i+j)
                if seq_read_R[i+j] != seq_next_read_R[i+start]:
                    localMismatches.add(i+j)
                # if len(mismatches) > m2:
                #     return False

        # Check if next_read_R matches pre-consensus
        for i in xrange(beginning):
            # print group.preConsensus
            if len(group.preConsensus[i+l]) > 1 or group.\
                   preConsensus[i+l].keys()[0] != seq_next_read_R[i+toExtend]:
                mismatches.add(beginning-(i+l))
                if seq_read_R[i+l] != seq_next_read_R[i+toExtend]:
                    localMismatches.add(beginning-(i+l))
                # if len(mismatches) > m2:
                #     return False

        if len(mismatches) > m2 and len(localMismatches) <= m2:
            createNewGroup(groupOrg, next_read_R, seq_next_read_R, offset, localMismatches)
            return False
        elif len(mismatches) > m2:
            return False

        #print mismatches
        group.mismatches = mismatches


    # Update pre-consensus
    for i in xrange(beginning):
        bp = seq_next_read_R[i+toExtend]
        group.preConsensus[i+l][bp] = group.preConsensus[i+l].get(bp, 0) + 1

    # Update main-consensus
    for i in xrange(len(seq_next_read_R)-start):
        bp = seq_next_read_R[i+start]
        group.consensus[i+j][bp] = group.consensus[i+j].get(bp, 0) + 1

    # Extend pre-consensus if required
    if toExtend > 0:
        prePart = [{seq_next_read_R[i]:1} for i in xrange(toExtend)]
        group.preConsensus = prePart + group.preConsensus

    # Add read to group
    global c2
    c2 += 1
    if next_read_R in group.rightParts:
        group.rightParts[next_read_R].add(offset)
    else:
        group.rightParts[next_read_R] = set([offset])

    #groupOrg.mismatches += mismatches
    return True


def dismissRead(group, seqs, next_read_R, offset):
    seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
    #print " "*-offset,
    # print_fullConsensus([],group.consensus)
    # print " "*offset, seq_next_read_R
    # print
    lenToCompare = len(seq_next_read_R) - (group.leftReadsOffset - offset)
    mismatches = set()
    for i in xrange(lenToCompare):
        if len(group.consensus[group.leftReadsOffset+i]) > 1 or \
                group.consensus[group.leftReadsOffset+i].keys()[0] != \
                seq_next_read_R[i+group.leftReadsOffset-offset]:
            mismatches.add(group.leftReadsOffset+i)
            # if len(mismatches) > m1:
            #     return True

    # if offset > 0:
    #     lenToCompare = group.leftReadsOffset - offset
    #     j = offset
    # else:
    #     lenToCompare = group.leftReadsOffset
    #     l = -offset
    #
    # for i in xrange(lenToCompare):


    return mismatches


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
    #     print next_read_R
    #     print "LOL-NEW"
    #     print group.mismatches
    #     print newGroup.mismatches
    #     print mismatches
    #     print

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
    newGroup.rightParts[next_read_R] = set([offset])
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
    #     if seq_read_R[i+j] != seq_next_read_R[i+start]:
    #         all_mismatches.add(i+j)
    #         if len(all_mismatches) > m2:
    #             return True

    # Check consensus
    # start = beginning + toExtend
    # for i in xrange(group.leftReadsOffset-j):
    #     if len(rightPartGroup.consensus[i+j]) > 1 or \
    #            rightPartGroup.consensus[i+j].iterkeys().next() \
    #            != seq_next_read_R[i+start]:
    #         all_mismatches.add(i+j)
    #         if len(all_mismatches) > m2:
    #             return False


    # Check pre-consensus
    for i in xrange(beginning):
        if len(rightPartGroup.preConsensus[i+l]) > 1 or \
               rightPartGroup.preConsensus[i+l].iterkeys().next() \
               != seq_next_read_R[i+toExtend]:
            all_mismatches.add((i+l)-len(rightPartGroup.preConsensus))
            if len(all_mismatches) > m2:
                # if next_read_R == 11035:
                #     print lenPreConsensus
                #     print "hej"
                #     print offset
                #     print toExtend
                #     print beginning
                #     print all_mismatches
                #     print_fullConsensus(rightPartGroup.preConsensus, rightPartGroup.consensus)
                #     print next_read_R
                # createNewGroup(group, next_read_R, seq_next_read_R, offset,
                #                mismatches)
                return False


    # Add existing mismatches
    # for mis in group.mismatches:
    #     all_mismatches.add(mis)
    for mis in rightPartGroup.mismatches:
        all_mismatches.add(mis)
        if len(all_mismatches) > m2:
            # createNewGroup(group, next_read_R, seq_next_read_R, offset,
            #                mismatches)
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
    #     print next_read_R
    #     print "LOL"
    #     print group.mismatches
    #     print rightPartGroup.mismatches
    #     print
    if next_read_R in rightPartGroup.rightParts:
        rightPartGroup.rightParts[next_read_R].add(offset)
    else:
        rightPartGroup.rightParts[next_read_R] = set([offset])

    rightPartGroup.mismatches = all_mismatches
    return True


def testRead(group, seqs, read_R, next_read_R, offset, m2, log):
    seq_read_R = seqs[read_R-1]+seqs[read_R]
    seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]

    # Check overlapping part
    lenToCompare = len(seq_next_read_R) - (group.leftReadsOffset - offset)
    mismatches = set([mis for mis in group.mismatches])
    for i in xrange(lenToCompare):
        # if next_read_R == 1805:
        #     print "hejhej"
        #     print group.consensus[group.leftReadsOffset+i].keys()[0], seq_next_read_R[i+group.leftReadsOffset-offset]
        #     print group.leftReadsOffset+i
        #     print group.leftReadsOffset
        if len(group.consensus[group.leftReadsOffset+i]) > 1 or \
                group.consensus[group.leftReadsOffset+i].iterkeys().next() \
                != seq_next_read_R[i+group.leftReadsOffset-offset]:
            mismatches.add(group.leftReadsOffset+i)
            # if next_read_R == 1805:
            #     #print " "*-offset,
            #     print_fullConsensus([], group.consensus)
            #     print " "*offset, seq_next_read_R
            #     print offset
            #     print "lololol"
            #     print group.consensus[group.leftReadsOffset+i].keys()[0], seq_next_read_R[i+group.leftReadsOffset-offset]
            if len(mismatches) > m2:
                return False

    # Check rest against anchor
    if offset > 0:
        j = offset
        k = 0
    else:
        j = 0
        k = -offset
    for i in xrange(group.leftReadsOffset-abs(offset)):
        if seq_read_R[i+j] != seq_next_read_R[i+k]:
            mismatches.add(i+j)
            if len(mismatches) > m2:
                return False

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


def alignRightParts(read_R, seqs, alignedGroups, candidatePairs, log):
    seen = set()
    for group in alignedGroups:
        ids = []
        for read_L in group.leftParts:
            for next_read_R in candidatePairs[read_L]:
                if read_R != next_read_R and \
                             next_read_R not in group.checkedRightParts:
                    if next_read_R == 2119:
                        print next_read_R
                        print_fullConsensus([], group.consensus)
                        print seqs[next_read_R-1]+seqs[next_read_R]
                    for alignInfo in findAlignment(seqs[next_read_R],
                                    seqs[read_L], group.readROffset, 0, next_read_R):
                        for offset2 in group.leftParts[read_L]:
                            ids.append(next_read_R)
                            offset1, mis, lenCompared = alignInfo
                            if next_read_R == 2119:
                                print alignInfo
                            global numreadR
                            numreadR += 1
                            offset = offset2 - offset1
                            offset += group.readROffset - \
                                      len(seqs[next_read_R-1])
                            
                            # print_fullConsensus([], group.consensus)
                            # print " "*-offset, seqs[read_R-1]+seqs[read_R]
                            # print " "*offset, seqs[next_read_R-1]+seqs[next_read_R]
                            #
                            # print "offset", offset

                            # mis = dismissRead(group, seqs, next_read_R,offset)
                            #
                            # if len(mis) > M2:
                            #     group.checkedRightParts.add(next_read_R)
                            #     continue
                            #
                            # added = False
                            # for rightPartGroup in group.rightPartGroups:
                            #     added = fitsInGroup4(group, rightPartGroup,
                            #         seqs, read_R, next_read_R, offset,
                            #         mis, M2)
                            #
                            # # if not added:
                            # if len(group.rightPartGroups) == 0:
                            #     seq_next_read_R = seqs[next_read_R-1]+seqs[next_read_R]
                            #     createNewGroup(group, next_read_R,
                            #             seq_next_read_R, offset, mis)

                            testRead(group, seqs, read_R, next_read_R,
                                     offset, M2, log)
                            group.checkedRightParts.add(next_read_R)
                            #ids.append(next_read_R)
                #group.checkedRightParts.add(next_read_R)
        #print sorted(ids)
        # ids2 = [(x,y) for x, y in Counter(ids).items() if y > 1]
        # if len(ids2) > 0 and read_R not in seen:
        #     print read_R,
        #     #print ids2
        #     seen.add(read_R)
        # print group.checkedRightParts
        # print ids



        # Second pass to add possible missing pre-right-parts to other groups
        # alignRightParts_secondPass(group, seqs)
        # for shallowGroup in group.shallowGroups:
        #     alignRightParts_secondPass(shallowGroup, seqs)

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

    # For testing only
    # input_file = "candidate_pairs_k_16_b_2_r_5_m_6.txt"

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

        return (time.clock() - totim)

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
    """
    Main method call
    """
    main()

    sys.exit(0)
