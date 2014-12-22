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
                      help="<VALUE> defines the minhash algorithm used.")

    parser.add_option("-s", "--seed",
                      metavar="<VALUE>",
                      type=int,
                      action="store",
                      dest="s",
                      help="set <VALUE> as seed for hash functions.")

    parser.add_option("-l", "--log_file",
                      metavar="<VALUE>",
                      type=str,
                      default="log.txt",
                      action="store",
                      dest="log",
                      help="set <VALUE> as seed for hash functions.")

    (options, args) = parser.parse_args()

    return options.fasta_file, options.k, options.threshold,\
        options.bands, options.rows, options.sim, options.m, \
        options.s, options.log


def memory_usage_resource():
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
    for element in output:
        log_output += str(element) + " "
    print log_output
    log_output += "\n"
    log_file.write(log_output)
    if flush:
        log_file.flush()

# ************************************************************************** #
#                                                                            #
#                              Pre-computations                              #
#                                                                            #
# ************************************************************************** #
def computeHashFunctions(n, shingles, log):
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
    Extract the reads (DNA sequences) from the given fasta file
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
        numPairs = 0

        # Minhashing using pre-computed hash functions
        if minhash_alg == 3 or minhash_alg > 5:
            shingles = computeShinglesTable(fasta_file, k, log)
        else:  # minhash alg 1, 2, 4 or 5
            shingles = computeShinglesSet(fasta_file, k, log)
        for b in xrange(bands):
            buckets = dict()
            minhashing(fasta_file, shingles, buckets, k, rows,
                       minhash_alg, b, bands, log)
            lshBand(buckets, b, candidatePairs, log)


        # Type 1 - "First Hash"
        # if minhash_alg < 3:
        #     shingles = computeShinglesSet(fasta_file, k, log)
        #     for b in xrange(bands):
        #         hashfuncs = computeHashFunctions(rows, shingles, log)
        #         buckets = dict()
        #         idx = 0
        #         for part in getPartsFromFile(fasta_file, log):
        #             minhashing_alg1(part, idx, shingles, buckets, k, rows,
        #                             hashfuncs)
        #             idx += 1
        #         lshBand(buckets, b, candidatePairs, log)


        logprint(log, False, "\nNumber of unique candidate pairs",
                 sum(len(candidatePairs[i]) for i in candidatePairs))
        logprint(log, False, "Finished LSH in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())

        return candidatePairs

    else:
        logprint(log, True, "ERROR: NO FASTA FILE PROVIDED")
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
        else: # Default to minhash alg 6
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

    logprint(log, True, "Number of buckets in band", b, ":", len(buckets))
    numPairs = 0
    for bucket in buckets:
        inBucket = buckets[bucket]
        numPairs += len(inBucket) * (len(inBucket)-1) / 2
    logprint(log, False, "Number of candidate pairs in band", b,
             ":", numPairs)
    logprint(log, True, "Number of unique candidate pairs in band", b,
             ":", numPairsUnique)

        # print "Finished LSH for band", b, "in", (time.clock() - tim) / 60, \
        #       "minutes"

    return candidatePairs


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
    enumerates them, then store them in table for fast look up.
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
def sequenceAlignment(candidatePairs, fasta_file, log):
    seqs = getAllReads(fasta_file, log)
    alignMatrix = dict()
    alignedLeftGroups = dict()
    alignedRightGroups = dict()

    numParts = len(candidatePairs) / 2
    proc = 0
    tim = time.clock()
    for read_R in candidatePairs:
        if read_R % 2 == 1:
        #if read_R == 275:

            # Align left parts
            alignLeftParts(read_R, seqs, alignedLeftGroups, alignMatrix,
                           candidatePairs, log)

            # Align right parts
            newOffset = alignRightParts(read_R, seqs, alignedLeftGroups,
                                        alignedRightGroups, alignMatrix,
                                        candidatePairs, log)

            logprint(log, False, "\n", " "*(abs(newOffset)), seqs[read_R-1]+" "+seqs[read_R])
            for group in alignedLeftGroups[read_R]:
                for read_L in group:
                    if read_L % 2 == 0:
                        logprint(log, False, " " * (abs(newOffset) + 2 +
                             alignMatrix[read_L] + len(seqs[read_R-1])),
                             seqs[read_L]+" "+seqs[read_L+1])
                    else:
                        logprint(log, False, " " * (abs(newOffset) +
                             alignMatrix[read_L] + 1),
                             seqs[read_L-1]+" "+seqs[read_L])
                logprint(log, False, sorted(list(group)))
                logprint(log, False, len(group), "\n")
            logprint(log, False, " "*(abs(newOffset)), seqs[read_R-1]+
                                 " "+seqs[read_R])
            for group in alignedRightGroups[read_R]:
                for read_R in group:
                    logprint(log, False, " " * (abs(newOffset) +
                             alignMatrix[read_R]), seqs[read_R-1]+
                             " "+seqs[read_R])
                logprint(log, False, sorted(list(group)))
                logprint(log, False, len(group))

            #sys.exit()

            proc += 1
            if proc % 500000 == 0:
                logprint(log, True, "Processed", proc, "of", numParts,
                         "right parts in", (time.clock()-tim) / 60, "minutes")

    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


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


def alignLeftParts(read_R, seqs, leftGroups, alignMatrix, candidatePairs, log):
    leftGroups[read_R] = []
    for read_L in candidatePairs[read_R]:
        offset = findAlignment(seqs[read_R], seqs[read_L], 0, log)
        if offset > -1:
            alignMatrix[read_L] = offset
            if len(leftGroups[read_R]) == 0:
                leftGroups[read_R].append(set([read_L]))
            else:
                newGroup = True
                for group in leftGroups[read_R]:
                    if fitsInGroup(group, seqs, read_R, read_L,
                                   alignMatrix, 0):
                        group.add(read_L)
                        newGroup = False
                        #break
                if newGroup:
                    leftGroups[read_R].append(set([read_L]))

    # if len(alignedGroups[read_R]) > 1:
    #     print "id", read_R, "cand. pairs", len(candidatePairs[read_R])
    #     for group in alignedGroups[read_R]:
    #         print group
    #         #print len(group)
    #
    #     for group in alignedGroups[read_R]:
    #         print seqs[read_R]
    #         for read_L in group:
    #             print " "*(alignMatrix[read_L]-1), seqs[read_L]
    #         print

def findAlignment(read_R, read_L, allowed_mismatches, log):
    # print "read_R:", read_R
    # print "read_L:", read_L
    # print
    offset = 0
    doPrint = False
    if len(read_R) > len(read_L):
        offset = len(read_R) - len(read_L)
    lengthToCompare = len(read_L)
    while lengthToCompare > 10:
        mismatches = 0
        for i in xrange(lengthToCompare-1):
            if read_R[i+offset] != read_L[i]:
                mismatches += 1
                if mismatches > allowed_mismatches:
                    break
        if mismatches <= allowed_mismatches:
            if doPrint:
                print "offset:", offset
                print read_R
                if offset == 0:
                    print read_L
                else:
                    print " "*(offset-1), read_L
            return offset
        offset += 1
        lengthToCompare -= 1
    return -1


def fitsInGroup(group, seqs, read_R, read_L, alignMatrix, al_mismatches):
    # if group == []:
    #     return False
    lread_R = seqs[read_R]
    lread_L = seqs[read_L]+seqs[read_L+1]
    #offset1 = alignmentMatrix[read_R][read_L]
    offset1 = alignMatrix[read_L]
    extraPart1 = len(lread_L)-len(lread_R)+offset1
    for read in group:
        #offset2 = alignmentMatrix[read_R][read]
        lread = seqs[read]+seqs[read+1]
        offset2 = alignMatrix[read]
        extraPart2 = len(lread)-len(lread_R)+offset2
        mismatches = 0
        for i in xrange(min(extraPart1, extraPart2)):
            if lread_L[i+len(lread_R)-offset1] != lread[i+len(lread_R)-offset2]:
                mismatches += 1
                if mismatches > al_mismatches:
                    return False
    return True


def alignRightParts(read_R, seqs, leftGroups, rightGroups, alignMatrix, candidatePairs, log):
    rightGroups[read_R] = []
    newOffset = 0
    for group in leftGroups[read_R]:
        #rightPartsLeft = set()
        #rightPartsRight = set()
        for read_L in group.copy():
            for next_read_R in candidatePairs[read_L]:
                if read_R != next_read_R and \
                             next_read_R not in alignMatrix:
                    offset = findAlignment(seqs[next_read_R],
                                           seqs[read_L], 0, log)
                    if offset > -1:
                        offset = alignMatrix[read_L] - offset
                        if offset < newOffset:
                            newOffset = offset
                        alignMatrix[next_read_R] = offset
                        if offset > 0:
                            if len(leftGroups[read_R]) == 0:
                                leftGroups[read_R].append(set([next_read_R]))
                            else:
                                newGroup = True
                                for group in leftGroups[read_R]:
                                    if fitsInLeftGroup(group, seqs,
                                       read_R, next_read_R, alignMatrix, 0):
                                           group.add(next_read_R)
                                           newGroup = False
                                           #break
                                if newGroup:
                                    leftGroups[read_R].\
                                    append(set([next_read_R]))
                        elif offset < 0:
                            if len(rightGroups[read_R]) == 0:
                                rightGroups[read_R].append(set([next_read_R]))
                            else:
                                newGroup = True
                                for group in rightGroups[read_R]:
                                    if fitsInRightGroup(group, seqs,
                                       read_R, next_read_R, alignMatrix, 0):
                                           group.add(next_read_R)
                                           newGroup = False
                                           #break
                                if newGroup:
                                    rightGroups[read_R].\
                                    append(set([next_read_R]))

                            # elif fitsInRightGroup(rightPartsLeft, seqs,
                            #         read_R, next_read_R, alignMatrix, 0):
                            #     rightGroup.add(next_read_R)
                            # else:
                            #     print logprint(log, False, " "*(abs(offset)),
                            #                    seqs[read_R])
                            #     print logprint(log, False, seqs[read_R])


                        # print alignMatrix[read_L]
                        # print offset
                        # if alignMatrix[read_L] > offset:
                        #     newOffset = alignMatrix[read_L] - offset
                        #     print seqs[read_R]
                        #     print " "*(newOffset-1), seqs[next_read_R]
                        #     print " "*(alignMatrix[read_L]-1), seqs[read_L]
                        # else:
                        #     newOffset = offset - alignMatrix[read_L]
                        #     print " "*(newOffset-1), seqs[read_R]
                        #     print seqs[next_read_R]
                        #     print " "*(offset-1), seqs[read_L]

    return newOffset

def fitsInLeftGroup(group, seqs, read_R, next_read_R, alignMatrix, m2):
    lread_R = seqs[read_R-1]+seqs[read_R]
    lread_L = seqs[next_read_R-1]+seqs[next_read_R]
    offset1 = alignMatrix[next_read_R]
    extraPart1 = len(lread_L)-len(lread_R)+offset1
    for read in group:
        if read % 2 == 0:
            lread = seqs[read]+seqs[read+1]
            offset2 = alignMatrix[read]+len(seqs[read_R-1])
        else:
            lread = seqs[read-1]+seqs[read]
            offset2 = alignMatrix[read]
        extraPart2 = len(lread)-len(lread_R)+offset2
        
        # print read_R, next_read_R, read
        # print read_R-1, next_read_R-1, read-1
        # print offset1, offset2
        # print seqs[read_R-1]+" "+seqs[read_R]
        # print " "*(offset1-1), seqs[next_read_R-1]+" "+seqs[next_read_R]
        # if read % 2 == 0:
        #     print " "*(offset2), seqs[read]+" "+seqs[read+1]
        # else:
        #     print " "*(offset2-1), seqs[read-1]+" "+seqs[read]
        #sys.exit()
        
        mismatches = 0
        for i in xrange(min(extraPart1, extraPart2)):
            if lread_L[i+len(lread_R)-offset1] != lread[i+len(lread_R)-offset2]:
                mismatches += 1
                if mismatches > m2:
                    print "lol"
                    return False
    return True


def fitsInRightGroup(group, seqs, read_R, next_read_R, alignMatrix, m2):
    for read in group:
        extraPart = min(abs(alignMatrix[read]), abs(alignMatrix[next_read_R]))
        offset = alignMatrix[read] - alignMatrix[next_read_R]

        offset1 = 0 # next_read_R
        offset2 = 0 # read
        if offset > 0:
            offset1 = offset
        elif offset < 0:
            offset2 = abs(offset)

        # print read_R, next_read_R, read
        # print read_R-1, next_read_R-1, read-1
        # print offset1, offset2
        # padding = max(abs(alignMatrix[read]), abs(alignMatrix[next_read_R]))
        # print " "*(padding), seqs[read_R-1]+" "+seqs[read_R]
        # print " "*(offset2), seqs[next_read_R-1]+" "+seqs[next_read_R]
        # print " "*(offset1), seqs[read-1]+" "+seqs[read]
        mismatches = 0
        seqRead = seqs[read-1]+seqs[read]
        seqNextRead = seqs[next_read_R-1]+seqs[next_read_R]
        for i in xrange(extraPart):
            if seqNextRead[i+offset1] != seqRead[i+offset2]:
                mismatches += 1
                if mismatches > m2:
                    return False
    return True


# GTCAA AGTTCAG
#          TCAGAA TGCCC
#         TTCAGAA TGCC

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
        seed, log_file = optionParse()

    # n (number of hash functions = length of minhash signature) is given by
    n = bands * rows

    with open(log_file, 'w') as log:

        candidatePairs = runLSH(fasta_file, bands, rows, n, k, seed,
                                minhash_alg, log)

        tim = time.clock()
        sequenceAlignment(candidatePairs, fasta_file, log)
        logprint(log, False, "Finished sequence alignment",
                 (time.clock() - tim) / 60, "minutes")

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
