#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from __future__ import print_function

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "11/09/2014"
__version__ = "$Revision: 1.1"

from optparse import OptionParser
from operator import itemgetter
from collections import Counter
import itertools
import sys
import time
import random
import SuperFastHash as sfh
import mmh3
# import os


def read_data_seq(fasta_file):
    with open(fasta_file, "rU") as fasta_file:
        read = ""
        tim = time.clock()
        print "Fasta file..."
        reads = []
        for line in fasta_file:
            # If line starts with ">", which indicates end of a sequence, append it to list of reads
            if line.startswith(">"):
                if read != "":
                    # Splits the string into two parts
                    leftpart = read[:len(read)/2]
                    rightpart = read[len(read)/2:]
                    reads.append(leftpart)
                    reads.append(rightpart)
                    read = ""
            # Concatenate multi-line sequences into one string
            else:
                read += line.strip().upper()

        print "Finished reading in", (time.clock() - tim) / 60, "minutes"
        print "Memory usage (in mb):", memory_usage_resource()


def LSH_run(fasta_file, bands, rows, n, k, seed, minhash_alg, log):
    if fasta_file:
        tim = time.clock()

        # random.seed(seed)
        # j = superFastHash(["hej","med","a","s","f"])
        # for i in j:
        #     print i
        # sys.exit()

        shingles = computeAllShingles(fasta_file, k, log)

        hashfuncs = None
        precomputeHash = True
        if precomputeHash:
            hashfuncs = computeHashes(n, shingles, log)
        band_buckets = [dict() for _ in xrange(bands)]
        read_minhash(fasta_file, shingles, band_buckets, k, rows, bands, n,
                           seed, hashfuncs, minhash_alg, log)
        LSH_band(band_buckets, log)

        # for b in xrange(bands):
        #     buckets = dict()
        #     idx = read_minhash(fasta_file, shingles, buckets, k, rows,
        #                        seed, hashfuncs)
        #     LSH_band(buckets, b)

        logprint(log, False, "Finished LSH in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())

    else:
        logprint(log, True, "ERROR: NO FASTA FILE GIVEN")
        sys.exit()


def computeHashes(n, shingles, log):
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    logprint(log, True, "Computed hashfunctions")
    return hashfuncs


def hashFuncGenerator2(shingles):
    h = range(len(shingles))
    random.shuffle(h)
    for i in h:
        yield i


def superFastHash(shingles):
    numShingles = len(shingles)
    for shingle in shingles:
        salt = random.randrange(numShingles)
        #print "salt:", salt
        yield sfh.SuperFastHash(shingle, salt) % numShingles


def superFastHashList(shingles, n):
    numShingles = len(shingles)
    for i in xrange(n):
        hashfunc = []
        for shingle in shingles:
            salt = random.randrange(numShingles)
            hashfunc.append(sfh.SuperFastHash(shingle, salt) % numShingles)
        yield hashfunc


def computeAllShingles(fasta_file, k, log):
    logprint(log, True, "Computing all shingles...")
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


def read_minhash(fasta_file, shingles, band_buckets, k, rows, bands, n, seed, hashfuncs, minhash_alg, log):
    with open(fasta_file, "rU") as fasta_file:
        read = ""
        tim = time.clock()
        logprint(log, True, "Minhashing...")
        idx = 0
        for line in fasta_file:
            # If line starts with ">", which indicates end of a sequence, append it to list of reads
            if line.startswith(">"):
                if read != "":
                    # Splits the string into two parts
                    leftpart = read[:len(read)/2]
                    if minhash_alg == 1:
                        minhashing_run1(leftpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs)
                    elif minhash_alg == 2:
                        minhashing_run2(leftpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs)
                    elif minhash_alg == 3:
                        shinglePos = {shingle: i for i, shingle in
                                      enumerate(shingles)}
                        minhashing_run3(leftpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs,
                                   shinglePos)
                    idx += 1

                    rightpart = read[len(read)/2:]
                    if minhash_alg == 1:
                        minhashing_run1(rightpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs)
                    if minhash_alg == 2:
                        minhashing_run2(rightpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs)
                    if minhash_alg == 3:
                        shinglePos = {shingle: i for i, shingle in
                                      enumerate(shingles)}
                        minhashing_run3(rightpart, idx, shingles, band_buckets,
                                   k, rows, bands, n, seed, hashfuncs,
                                   shinglePos)
                    idx += 1
                    read = ""

                if idx % 1000 == 0 and idx != 0:
                    logprint(log, True, "Processed", idx, "documents in",
                          (time.clock() - tim) / 60, "minutes")
            # Concatenate multi-line sequences into one string
            else:
                read += line.strip().upper()
        # Append the last read in the file to the list
        if read != "":
            leftpart = read[:len(read)/2]
            if minhash_alg == 1:
                minhashing_run1(leftpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs)
            elif minhash_alg == 2:
                minhashing_run2(leftpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs)
            elif minhash_alg == 3:
                shinglePos = {shingle: i for i, shingle in
                              enumerate(shingles)}
                minhashing_run3(leftpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs,
                           shinglePos)

            idx += 1
            rightpart = read[len(read)/2:]
            if minhash_alg == 1:
                minhashing_run1(rightpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs)
            if minhash_alg == 2:
                minhashing_run2(rightpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs)
            if minhash_alg == 3:
                shinglePos = {shingle: i for i, shingle in
                              enumerate(shingles)}
                minhashing_run3(rightpart, idx, shingles, band_buckets,
                           k, rows, bands, n, seed, hashfuncs,
                           shinglePos)
            idx += 1

            if idx % 1000 == 0:
                logprint(log, True, "Processed", idx, "documents in",
                      (time.clock() - tim) / 60, "minutes")

        logprint(log, False, "Finished minhashing in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def minhashing_run1(dna, idx, shingles, band_buckets, k, rows, bands, n, seed, hashfuncs):
    """
    Create minhash signatures using the shingles
    """
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    random.seed(seed)
    docShingles = getDocShingles(dna, k)
    signature = [None for i in xrange(n)]
    # For each row in the 'character matrix'
    for sigPos in xrange(n):
        for i, h in enumerate(hashfuncs[sigPos]):
        #for i, h in enumerate(hashFuncGenerator2(shingles)):
        #for i, h in enumerate(superFastHash(shingles)):
        # for i in xrange(len(shingles)):
        #     h = random.randrange(len(shingles))
            #print h
            if shingles[h] in docShingles:
                signature[sigPos] = i
                break
        if signature[sigPos] == None:
            signature[sigPos] = len(shingles)

    i = 0
    for buckets in band_buckets:
        # print signature[i:i+rows]
        key = int(''.join(map(str, signature[i:i+rows])))
        if key in buckets:
            buckets[key].append(idx)
        else:
            buckets[key] = [idx]
        i += rows


def minhashing_run2(dna, idx, shingles, band_buckets, k, rows, bands, n, seed, hashfuncs):
    """
    Create minhash signatures using the shingles
    """
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    random.seed(seed)
    docShingles = getDocShingles(dna, k)
    signature = [None for i in xrange(n)]
    # For each row in the 'character matrix'
    for r in xrange(len(shingles)):
        # If the shingle is in the document, then
        # if doc.id in shinglesDict[shingles[r]]:
        if shingles[r] in docShingles:
            # Find the 'first' shingle relative to each permutation
            # for i in xrange(n):
            #     if signature[i] is None or signature[i] > hashfuncs[i][r]:
            #         signature[i] = hashfuncs[i][r]
            # for i, hashfunc in enumerate(hashFuncGenerator(n, shingles)):
            #     if signature[i] is None or signature[i] > hashfunc[r]:
            #         signature[i] = hashfunc[r]
            for i, hashfunc in enumerate(superFastHashList(shingles, n)):
                if signature[i] is None or signature[i] > hashfunc[r]:
                    signature[i] = hashfunc[r]

    i = 0
    for buckets in band_buckets:
        key = int(''.join(map(str, signature[i:i+rows])))
        if key in buckets:
            buckets[key].append(idx)
        else:
            buckets[key] = [idx]
        i += rows


def minhashing_run3(dna, idx, shingles, band_buckets, k, rows, bands, n, seed, hashfuncs, shinglePos):
    """
    Create minhash signatures using the shingles
    """
    numShingles = len(shingles)

    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    random.seed(seed)
    signature = []
    docShingles = getDocShingles(dna, k)
    #for h in hashfuncs:
    #for h in hashFuncGenerator(n, shingles):
    for h in superFastHashList(shingles, n):
        minVal = numShingles+1
        for shingle in docShingles:
            pos = shinglePos[shingle]
            if h[pos] < minVal:
                minVal = h[pos]
        signature.append(minVal)
    # print signature
    # if None in signature:
    #     sys.exit()

    i = 0
    for buckets in band_buckets:
        key = int(''.join(map(str, signature[i:i+rows])))
        if key in buckets:
            buckets[key].append(idx)
        else:
            buckets[key] = [idx]
        i += rows


def LSH_band(band_buckets, log):
    tim = time.clock()
    logprint(log, True, "Running LSH...")
    b = 0
    candidatePairs = set()
    for buckets in band_buckets:
        numPairsUnique = 0
        b += 1
        for bucket in buckets:
            for i in xrange(len(buckets[bucket])):
                id1 = buckets[bucket][i]
                for j in xrange(i+1, len(buckets[bucket])):
                    id2 = buckets[bucket][j]
                    if id1 % 2 == 0 and id2 % 2 == 1:
                        if id1 != id2 + 1:
                            candidatePairs.add((id1,id2))
                            numPairsUnique += 1
                    if id1 % 2 == 1 and id2 % 2 == 0:
                        if id1 - 1 != id2:
                            candidatePairs.add((id1,id2))
                            numPairsUnique += 1

        logprint(log, True, "Number of buckets in band", b, ":", len(buckets))
        numPairs = 0
        for bucket in buckets:
            inBucket = buckets[bucket]
            numPairs += len(inBucket) * (len(inBucket)-1) / 2
        logprint(log, False, "Number of candidate pairs in band", b,
                 ":", numPairs)
        logprint(log, True, "Number of candidate pairs in band", b,
                 ":", numPairsUnique)

        # print "Finished LSH for band", b, "in", (time.clock() - tim) / 60, \
        #       "minutes"

    logprint(log, False, "Number of unique candidate pairs",
             len(candidatePairs))


# *************************************************************************** #
#                                                                             #
#                                Data structure                               #
#                                                                             #
# *************************************************************************** #
class Document(object):
    """
    Defines the data structure of a document. Like dna sequence, id, left or
    right part of string, shingles related to the string, the minhash signature
    """
    dna = ""
    id = 0
    isLeft = False  # If document comes from the left or right part of a read
    # shingles = set()
    # counterShingles = None
    signature = []
    similarTo = set()

    # Initializer
    def __init__(self, dna, id, isLeft):
        self.dna = dna
        self.id = id
        self.isLeft = isLeft
        self.similarTo = set()


def createDocument(dna, id, isLeft):
    """
    Function to create document objects
    """
    document = Document(dna, id, isLeft)
    return document


# *************************************************************************** #
#                                                                             #
#                                   Helpers                                   #
#                                                                             #
# *************************************************************************** #
def optionParse():
    """
    Parse arguments from command line.
    """
    desc = """Compare sets of sequencing data to find mutations."""

    parser = OptionParser(usage="usage: %prog --fasta_file filename",
                          description=desc,
                          version="%prog version 1.0")

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
                      default="1",
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


def readDataChunks(fasta_file, log):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    if fasta_file:
        logprint(log, False, "Processing fasta file...")
        tim = time.clock()
        documents = []
        seqs = 0
        id = 0
        with open(fasta_file, "rU") as fasta_file:
            read = ""
            for line in fasta_file:
                # If line starts with ">", which indicates end of a sequence
                # Then, append it to list of reads
                if line.startswith(">"):
                    if read != "":
                        # Splits the string into two parts
                        leftpart = read[:len(read)/2]
                        rightpart = read[len(read)/2:]
                        # Creates a doc object for the left part of the read
                        #documents.append(createDocument(leftpart, id, True))
                        yield createDocument(leftpart, id, True)
                        id += 1
                        # Creates a doc object for the right part of the read
                        #documents.append(createDocument(rightpart, id, False))
                        yield createDocument(rightpart, id, False)
                        id += 1
                        read = ""
                        seqs += 1
                # Concatenate multi-line sequences into one string
                else:
                    read += line.strip().upper()

                # if id != 0 and id % 1000 == 0:
                #     logprint(log, False, "Read", format(seqs, ',d'),
                #         "sequences in", (time.clock() - tim) / 60, "minutes")
                #     logprint(log, True, "Memory usage (in mb):",
                #              memory_usage_resource())
                #     yield documents
                #     seqs = 0
                #     documents = []

            # Append the last read in the file to the list
            if read != "":
                # Splits the string into two parts
                leftpart = read[:len(read)/2]
                rightpart = read[len(read)/2:]
                # Creates a document object for the left part of the read
                documents.append(createDocument(leftpart, id, True))
                id += 1
                # Creates a document object for the right part of the read
                documents.append(createDocument(rightpart, id, False))
                seqs += 1

        logprint(log, False, "Read", format(seqs, ',d'), "sequences in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
        yield documents

    else:
        logprint(log, False, "ERROR: NO FASTA FILE GIVEN")
        sys.exit()


def readData(fasta_file, log):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    if fasta_file:
        logprint(log, False, "Processing fasta file...")
        tim = time.clock()
        documents = []
        seqs = 0
        id = 0
        with open(fasta_file, "rU") as fasta_file:
            read = ""
            for line in fasta_file:
                # If line starts with ">", which indicates end of a sequence
                # Then, append it to list of reads
                if line.startswith(">"):
                    if read != "":
                        # Splits the string into two parts
                        leftpart = read[:len(read)/2]
                        rightpart = read[len(read)/2:]
                        # Creates a doc object for the left part of the read
                        documents.append(createDocument(leftpart, id, True))
                        id += 1
                        # Creates a doc object for the right part of the read
                        documents.append(createDocument(rightpart, id, False))
                        id += 1
                        read = ""
                        seqs += 1
                # Concatenate multi-line sequences into one string
                else:
                    read += line.strip().upper()
            # Append the last read in the file to the list
            if read != "":
                # Splits the string into two parts
                leftpart = read[:len(read)/2]
                rightpart = read[len(read)/2:]
                # Creates a document object for the left part of the read
                documents.append(createDocument(leftpart, id, True))
                id += 1
                # Creates a document object for the right part of the read
                documents.append(createDocument(rightpart, id, False))
                seqs += 1

        logprint(log, False, "Read", format(seqs, ',d'), "sequences in",
                 (time.clock() - tim) / 60, "minutes")
        logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
        return documents

    else:
        logprint(log, False, "ERROR: NO FASTA FILE GIVEN")
        sys.exit()


# *************************************************************************** #
#                                                                             #
#                          Locality Sensitive Hashing                         #
#                                                                             #
# *************************************************************************** #
def getDocShingles(dna, k):
    # shingles = set()
    # for i in xrange(len(doc.dna)-k+1):
    #     # create k-shingle (substring) from the document
    #     shingle = doc.dna[i:i+k]
    #     # Add it to the set of all shingles
    #     shingles.add(shingle)
    #shingles = set(doc.dna[i:i+k] for i in xrange(len(doc.dna)-k+1))
    shingles = {dna[i:i+k] for i in xrange(len(dna)-k+1)}
    return shingles


def shingling(docs, k, log):
    tim = time.clock()
    logprint(log, False, "Shingling...")

    shingles = {doc.dna[i:i+k] for doc in docs for i in
                xrange(len(doc.dna)-k+1)}
    # shingles = set()  # Contains all k-shingles in the dataset
    # for doc in docs:
    #     for i in xrange(len(doc.dna)-k+1):
    #         # create k-shingle (substring) from the document
    #         shingle = doc.dna[i:i+k]
    #         # Add it to the set of all shingles
    #         shingles.add(shingle)

    logprint(log, False, "Finished shingling in", (time.clock() - tim) / 60,
             "minutes")
    logprint(log, False, "Number of shingles:", len(shingles))
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())
    return list(shingles)


def hashFuncGenerator(n, shingles):
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        yield h


def minhashingNew(docs, shingles, n, k, seed, log):
    """
    Create minhash signatures using the shingles
    """
    tim = time.clock()
    print "Computing hashfunctions..."
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    random.seed(seed)
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    print "Computed hashfunctions in", (time.clock() - tim) / 60, "minutes"

    tim = time.clock()
    logprint(log, False, "Minhashing...")
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    count = 0
    for doc in docs:
        count += 1
        docShingles = getDocShingles(doc.dna, k)
        signature = [None for i in xrange(n)]
        # For each row in the 'character matrix'
        for sigPos in xrange(n):
            for i, h in enumerate(hashfuncs[sigPos]):
                if shingles[h] in docShingles:
                    signature[sigPos] = i
                    break
        doc.signature = signature
        # print signature
        if count % 1000 == 0:
            logprint(log, True, "Processed", count, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    logprint(log, False, "Finished minhashing in", (time.clock() - tim) / 60,
             "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def minhashingOld(docs, shingles, n, k, seed, log):
    """
    Create minhash signatures using the shingles
    """
    # Create n different permutations (hash functions) of the shingles
    random.seed(seed)
    hashfuncs = []
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    print "Computed hashfunctions"

    tim = time.clock()
    logprint(log, False, "Minhashing...")
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    count = 0
    for doc in docs:
        count += 1
        docShingles = getDocShingles(doc.dna, k)
        signature = [None for i in xrange(n)]
        # For each row in the 'character matrix'
        for r in xrange(len(shingles)):
            # If the shingle is in the document, then
            # if doc.id in shinglesDict[shingles[r]]:
            if shingles[r] in docShingles:
                # Find the 'first' shingle relative to each permutation
                # for i, hashfunc in enumerate(hashFuncGenerator(n, shingles)):
                #     if signature[i] is None or signature[i] > hashfunc[r]:
                #         signature[i] = hashfunc[r]
                for i in xrange(n):
                    if signature[i] is None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        #print signature
        doc.signature = signature
        if count % 1000 == 0:
            logprint(log, True, "Processed", count, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    logprint(log, False, "Finished minhashing in", (time.clock() - tim) / 60,
             "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def minhashingOld2(docs, shingles, n, k, seed, log):
    """
    Create minhash signatures using the shingles
    """
    tim = time.clock()
    # shinglePos = dict()
    # for i, shingle in enumerate(shingles):
    #     shinglePos[shingle] = i
    shinglePos = {shingle: i for i, shingle in enumerate(shingles)}
    print "dict time:", time.clock() - tim

    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    random.seed(seed)
    numShingles = len(shingles)
    for i in xrange(n):
        h = range(numShingles)
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    print "Computed hashfunctions"
    log.write("Computed hashfunctions\n")

    tim = time.clock()
    logprint(log, False, "Minhashing...")
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    #docShingleTime = 0
    count = 0
    for doc in docs:
        count += 1
        signature = []
        #tim2 = time.clock()
        docShingles = getDocShingles(doc.dna, k)
        #docShingleTime += time.clock() - tim2
        for h in hashfuncs:
            #minVal = min((h[shinglePos[shingle]] for shingle in docShingles))
            minVal = numShingles+1
            for shingle in docShingles:
                pos = shinglePos[shingle]
                if h[pos] < minVal:
                    minVal = h[pos]
            signature.append(minVal)
        doc.signature = signature
        # print signature
        if count % 1000 == 0:
            logprint(log, True, "Processed", count, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    #print "docShingleTime:", docShingleTime
    logprint(log, False, "Finished minhashing in", (time.clock() - tim) / 60,
             "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def minhashingOld3(docs, shingles, n, k, seed, log):
    """
    Create minhash signatures using the shingles
    """
    random.seed(seed)
    numShingles = len(shingles)
    hashfuncs = []
    for i in xrange(n):
        h = range(numShingles)
        random.shuffle(h)
        shinglePos = {shingle: h[i] for i, shingle in enumerate(shingles)}
        hashfuncs.append(shinglePos)

    tim = time.clock()
    logprint(log, False, "Minhashing...")
    count = 0
    for doc in docs:
        count += 1
        signature = []
        docShingles = getDocShinglesOld(doc, k)
        for h in hashfuncs:
            #print [h[shingle] for shingle in docShingles]
            #print min(h[shingle] for shingle in docShingles)
            #signature.append(min(h[shingle] for shingle in docShingles))
            minVal = numShingles+1
            for shingle in docShingles:
                if h[shingle] < minVal:
                    minVal = h[shingle]
            signature.append(minVal)
        doc.signature = signature
        if count % 1000 == 0:
            logprint(log, True, "Processed", count, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    logprint(log, False, "Finished minhashing in", (time.clock() - tim) / 60,
             "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def minhashing_mem_ef(docs, shingles, buckets, k, rows, seed, log):
    """
    Create minhash signatures using the shingles
    """
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    for i in xrange(rows):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    print "Computed hashfunctions"

    tim = time.clock()
    logprint(log, False, "Minhashing...")
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    count = 0
    for doc in docs:
        count += 1
        docShingles = getDocShingles(doc.dna, k)
        signature = [None for i in xrange(rows)]
        # For each row in the 'character matrix'
        for r in xrange(len(shingles)):
            # If the shingle is in the document, then
            # if doc.id in shinglesDict[shingles[r]]:
            if shingles[r] in docShingles:
                # Find the 'first' shingle relative to each permutation
                for i in xrange(rows):
                    if signature[i] is None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        # print sum(signature)
        key = int(''.join(map(str, signature)))
        if key in buckets:
            buckets[key].append(doc)
        else:
            buckets[key] = [doc]
        # print signature
        # if count < 20:
        #     print signature
        if count % 1000 == 0:
            logprint(log, True, "Processed", count, "documents in",
                     (time.clock() - tim) / 60, "minutes")

    logprint(log, False, "Finished minhashing in", (time.clock() - tim) / 60,
                 "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def LSH(documents, bands, rows, shingles, k, seed, log):
    """
    Use Locality-Sensitive Hashing (LSH) with the banding technique to
    hash similar elements to the same buckets.
    """
    index = 0
    random.seed(seed)
    for b in xrange(bands):
        # For each band, create a bucket array with signature as key
        buckets = dict([])
        minhashing_mem_ef(documents, shingles, buckets, k, rows, seed, log)

        tim = time.clock()
        logprint(log, False, "Running LSH...")
        numPairs = 0
        for bucket in buckets:
            for i in xrange(len(buckets[bucket])):
                doc1 = buckets[bucket][i]
                for j in xrange(i+1, len(buckets[bucket])):
                    doc2 = buckets[bucket][j]
                    if doc1.isLeft and not doc2.isLeft:
                        if doc1.id + 1 != doc2.id:
                            doc1.similarTo.add(doc2)
                            numPairs += 1
                    if not doc1.isLeft and doc2.isLeft:
                        if doc1.id != doc2.id + 1:
                            doc1.similarTo.add(doc2)
                            numPairs += 1

        logprint(log, False, "Number of buckets in band", b+1, ":",
                 len(buckets))
        # numPairs = 0
        # for bucket in buckets:
        #     inBucket = buckets[bucket]
        #     numPairs += len(inBucket) * (len(inBucket)-1) / 2
        logprint(log, False, "Number of candidate pairs in band", b+1, ":",
                 numPairs)

        logprint(log, True, "Finished LSH for band", b, "in",
                 (time.clock() - tim) / 60, "minutes")

    count = 0
    for doc in documents:
        count += len(doc.similarTo)
    logprint(log, False, "Number of unique candidate pairs", count)

    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


def LSH_old(documents, bands, rows, log):
    """
    Use Locality-Sensitive Hashing (LSH) with the banding technique to
    hash similar elements to the same buckets.
    """
    tim = time.clock()
    logprint(log, False, "Running LSH...")
    index = 0
    for b in xrange(bands):
        # For each band, create a bucket array with signature as key
        buckets = dict([])
        for doc in documents:
            # Obtain sub-signature of length rows
            col = doc.signature[index:index+rows]
            # Convert to string
            key = int(''.join(map(str, col)))
            # Place the document in a bucket
            if key in buckets:
                buckets[key].append(doc)
            else:
                buckets[key] = [doc]
        index += rows

        numPairsUnique = 0
        for bucket in buckets:
            for i in xrange(len(buckets[bucket])):
                doc1 = buckets[bucket][i]
                for j in xrange(i+1, len(buckets[bucket])):
                    doc2 = buckets[bucket][j]
                    #if doc1.isLeft and not doc2.isLeft:
                    if doc1.id % 2 == 0 and doc2.id % 2 == 1:
                        if doc1.id + 1 != doc2.id:
                            doc1.similarTo.add(doc2)
                            numPairsUnique += 1
                    #if not doc1.isLeft and doc2.isLeft:
                    if doc1.id % 2 == 1 and doc2.id % 2 == 0:
                        if doc1.id != doc2.id + 1:
                            doc1.similarTo.add(doc2)
                            numPairsUnique += 1
                    # if doc1.isLeft != doc2.isLeft:
                    #     doc1.similarTo.add(doc2)

        logprint(log, False, "Number of buckets in band", b+1, ":",
                 len(buckets))
        numPairs = 0
        for bucket in buckets:
            inBucket = buckets[bucket]
            numPairs += len(inBucket) * (len(inBucket)-1) / 2
        logprint(log, False, "Number of candidate pairs in band", b+1, ":",
                 numPairs)
        logprint(log, False, "Number of candidate pairs in band", b+1, ":",
                 numPairsUnique)

    count = 0
    for doc in documents:
        count += len(doc.similarTo)
    logprint(log, False, "Number of unique candidate pairs", count)
    logprint(log, False, "Finished LSH in", (time.clock()-tim) / 60, "minutes")
    logprint(log, True, "Memory usage (in mb):", memory_usage_resource())


# *************************************************************************** #
#                                                                             #
#                             Similarity checkers                             #
#                                                                             #
# *************************************************************************** #
def findSimilarPairs(documents, t, k, totim, numReads, sim_measure):
    """
    Find candidate pairs that has a similarity above the threshold t
    """
    counter = 0
    timer = time.clock()
    bestMatches1 = []
    bestMatches2 = []
    for doc1 in documents:
        for doc2 in doc1.similarTo:
            counter += 1
            if doc1.isLeft:
                bestMatch1 = globalAlignment(doc1, doc2, True)
                bestMatch2 = jaccardSim(doc1, doc2, k)
            else:
                bestMatch1 = globalAlignment(doc2, doc1, True)
                bestMatch2 = jaccardSim(doc2, doc1, k)
            bestMatches1.append(bestMatch1)
            bestMatches2.append(bestMatch2)

            if counter % 500000 == 0:
                print "Processed", format(counter2, ',d'), "pairs in", \
                      "time:", time.clock() - timer

    processing_time = time.clock() - timer
    c = "{:,}".format(counter).replace(",", ".")
    print "Processed", c, "pairs in time:", processing_time

    maxNumPairs = numReads * (numReads-1)
    with open("naive_vs_jaccard_standard_NA19240_2.txt", 'w') as f:
        f.write(str(numReads) + " " + str(maxNumPairs) + " " +
                str(counter) + " " + str(processing_time) + "\n")
        for idx, match1 in enumerate(bestMatches1):
            f.write(str(match1[0]) + " " + str(bestMatches2[idx]) + "\n")
            if idx % 500000 == 0:
                f.flush()

    bestMatches1.sort(key=itemgetter(0), reverse=False)
    with open("naive_NA19240_2.txt", 'w') as f:
        f.write(str(numReads) + " " + str(maxNumPairs) + " " +
                str(counter) + " " + str(processing_time) + "\n")
        for idx, match1 in enumerate(bestMatches1):
            f.write(str(idx) + " " + str(match1[0]) + " " +
                    str(match1[1]) + " " + str(match1[2]) + "\n")

    bestMatches2.sort()
    with open("standard_jaccard_NA19240_2.txt", 'w') as f:
        f.write(str(numReads) + " " + str(maxNumPairs) + " " +
                str(counter) + " " + str(processing_time) + "\n")
        for idx, match in enumerate(bestMatches2):
            f.write(str(idx) + " " + str(match) + "\n")


def findSimilarPairsOld(band_buckets, t, k, totim, output, numReads, sim_measure, writeToFile=False):
    """
    Find candidate pairs that has a similarity above the threshold t
    """
    counter = 0
    counter2 = 0
    timer = time.clock()
    bestMatches1 = []
    bestMatches2 = []
    for buckets in band_buckets:
        for bucket in buckets:
            if len(buckets[bucket]) > 1:
                for i in xrange(len(buckets[bucket])):
                    for j in xrange(i+1, len(buckets[bucket])):
                        doc1, doc2 = buckets[bucket][i], buckets[bucket][j]
                        counter += 1
                        # Check if doc1 and doc2 are candidate pairs
                        if doc1.isLeft != doc2.isLeft and doc1.id != doc2.id:
                            counter2 += 1
                            if counter2 < 10000000:
                                # euclideanDistance(doc1, doc2)
                                # testLCS(doc1, doc2)
                                # NeedlemanWunsch(doc1, doc2)
                                if doc1.isLeft:
                                    bestMatch1 = globalAlignment(doc1, doc2, True)
                                    bestMatch2 = jaccardSim(doc1, doc2, k)
                                    # if sim_measure == "naive":
                                    #     bestMatch = globalAlignment(doc1, doc2,
                                    #                                 True)
                                    # elif sim_measure == "jaccard":
                                    #     bestMatch = jaccardSim(doc1, doc2)
                                else:
                                    bestMatch1 = globalAlignment(doc2, doc1, True)
                                    bestMatch2 = jaccardSim(doc2, doc1, k)
                                    # if sim_measure == "naive":
                                    #     bestMatch = globalAlignment(doc2, doc1,
                                    #                                 True)
                                    # elif sim_measure == "jaccard":
                                    #     bestMatch = jaccardSim(doc2, doc1)

                                # bestMatches.append((bestMatch1, bestMatch2,
                                #                     doc1, doc2))
                                bestMatches1.append(bestMatch1)
                                bestMatches2.append(bestMatch2)

                                if counter2 % 500000 == 0:
                                    print "Processed", format(counter2, ',d') \
                                        .replace(",", "."), "pairs in", \
                                        "time:", time.clock() - timer
                            # else:
                            #     print "Total time used (in secs):", \
                            #         time.clock() - totim
                            #     bestMatches.sort(key=itemgetter(0),
                            #                      reverse=False)
                            #     with open("output.txt", 'w') as f:
                            #         counter = 0
                            #         for match in bestMatches:
                            #             f.write(str(counter) + "," +
                            #                     str(match[0])+"\n")
                            #             counter += 1
                            #     sys.exit()

    processing_time = time.clock() - timer
    c = "{:,}".format(counter2).replace(",", ".")
    print "Processed", c, "pairs in time:", processing_time
    print "{:,}".format(counter).replace(",", ".")
    print c
    #
    # maxNumPairs = numReads * (numReads-1)
    # with open("naive_vs_jaccard_standard_NA19240.txt", 'w') as f:
    #     f.write(str(numReads) + " " + str(maxNumPairs) + " " +
    #             str(counter2) + " " + str(processing_time) + "\n")
    #     for idx, match1 in enumerate(bestMatches1):
    #         f.write(str(match1[0]) + " " + str(bestMatches2[idx]) + "\n")
    #
    # bestMatches1.sort(key=itemgetter(0), reverse=False)
    # with open("naive_NA19240.txt", 'w') as f:
    #     f.write(str(numReads) + " " + str(maxNumPairs) + " " +
    #             str(counter2) + " " + str(processing_time) + "\n")
    #     for idx, match1 in enumerate(bestMatches1):
    #         f.write(str(idx) + " " + str(match1[0]) + " " +
    #                 str(match1[1]) + " " + str(match1[2]) + "\n")
    #
    # bestMatches2.sort()
    # with open("standard_jaccard_NA19240.txt", 'w') as f:
    #     f.write(str(numReads) + " " + str(maxNumPairs) + " " +
    #             str(counter2) + " " + str(processing_time) + "\n")
    #     for idx, match in enumerate(bestMatches2):
    #         f.write(str(idx) + " " + str(match) + "\n")

    # with open("output.txt", 'w') as f:
    #     for match in bestMatches:
    #         if match[0] > match[1] and match[0] != 1.0:
    #             f.write(str(match[0])+" "+str(match[1])+"\n")
    #             f.write(str(match[2].dna)+"\n")
    #             f.write(str(match[3].dna)+"\n")

    if writeToFile:
        if sim_measure == "naive":
            bestMatches.sort(key=itemgetter(0), reverse=False)
        elif sim_measure == "jaccard":
            bestMatches.sort()

        with open(output, 'w') as f:
            maxNumPairs = numReads*(numReads-1)/2
            f.write(str(numReads) + " " + str(maxNumPairs) + " " +
                    str(counter2) + " " + str(processing_time) + "\n")
            counter = 0
            for match in bestMatches:
                if sim_measure == "naive":
                    f.write(str(counter) + " " + str(match[0]) + " " +
                            str(match[1]) + " " + str(match[2]) + "\n")
                elif sim_measure == "jaccard":
                    f.write(str(counter) + " " + str(match) + "\n")
                counter += 1


def globalAlignment(doc1, doc2, extraInfo=False):
    """
    Aligning sequences by using a sliding window approach.
    Returns the best score (matches / seqlength) between the two sequences.
    """
    # dna1, dna2 = doc1.dna, doc2.dna
    # print doc1.dna
    # print len(doc1.dna)
    # print doc2.dna
    # print len(doc2.dna)
    # print

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

    # ## bp level jaccard sim - bag ## #
    # start = 0
    # bestScore = 0
    # seqLength = len(doc1.dna)-start
    # while seqLength > 10 and bestScore != 1.0:
    #
    #     # print seqLength
    #     # print doc1.dna
    #     # print doc2.dna
    #     # print doc1.dna[:seqLength]
    #     # print doc2.dna[-seqLength:]
    #
    #     # Count similar elements in signature
    #     counterA = Counter(doc1.dna[:seqLength])
    #     counterB = Counter(doc2.dna[-seqLength:])
    #     # Find intersection of elements in doc1 and doc2
    #     intersection = sum((counterA & counterB).values())
    #     # Find union of elements in doc1 and doc2
    #     #union = len(doc1.dna) + len(doc2.dna) - intersection
    #     union = 2 * seqLength - intersection
    #
    #     score = float(intersection) / union
    #     if score > bestScore:
    #         bestScore = score
    #
    #     start += 1
    #     seqLength = len(doc1.dna)-start
    #
    # return bestScore


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


# *************************************************************************** #
#                                                                             #
#                                Miscellaneous                                #
#                                                                             #
# *************************************************************************** #
def bucketSize(band_buckets):
    count = 0
    for buckets in band_buckets:
        for bucket in buckets:
            print len(buckets[bucket]),
            count += len(buckets[bucket])*(len(buckets[bucket])-1)/2
        print
        print count
        print
        # 712967 + 701401 + 690372 + 761370 + 686222 = 3552332
        # 3552332
        # 3916013
        # 3621369
        # 3804308

        # 3702466
        # 3742476
        # 3640087
        # 3768266


def compareAllPairs(documents, k, sim_measure):
    numReads = format(len(documents)/2, ',d').replace(",", ".")
    print "Number of reads:", numReads
    print "Number of docs:", format(len(documents), ',d').replace(",", ".")
    numPairs = len(documents)/2 * (len(documents)/2-1)
    strNumPairs = format(numPairs, ',d').replace(",", ".")
    print "Number of pairs to process:", strNumPairs
    bestMatches = [[] for i in xrange(37)]
    counter = 0
    timer = time.clock()
    for i in xrange(0, len(documents)):
        for j in xrange(i+1, len(documents)):
            doc1, doc2 = documents[i], documents[j]
            # Check if doc1 and doc2 are candidate pairs
            if doc1.isLeft != doc2.isLeft and doc1.id != doc2.id:
                counter += 1
                if doc1.isLeft:
                    bestMatch_naive = globalAlignment(doc1, doc2, True)
                    bestMatch_jaccard = jaccardSim(doc1, doc2, sim_measure)
                else:
                    bestMatch_naive = globalAlignment(doc2, doc1, True)
                    bestMatch_jaccard = jaccardSim(doc2, doc1, sim_measure)

                bestMatches[bestMatch_naive[2]-1].append((bestMatch_naive[0],
                                                          bestMatch_jaccard))
                # bestMatches.append(bestMatch_jaccard)

                # if bestMatch_naive[2] > 35 and bestMatch_jaccard > 0.4 and \
                #    bestMatch_jaccard < 0.8 and bestMatch_naive[0] > 0.8:
                #     print bestMatch_naive[0], bestMatch_jaccard
                #     print bestMatch_naive[1]
                #     print bestMatch_naive[2]
                #     print doc1.dna
                #     print doc2.dna
                #     print doc1.shingles.intersection(doc2.shingles)
                #     print doc1.shingles.union(doc2.shingles)

                if counter % 200000 == 0:
                    print "Processed", format(counter, ',d').replace(",", "."),\
                          "of", strNumPairs, "pairs in time (minutes):", \
                          (time.clock() - timer) / 60

    # with open("all_pairs_naive.txt", 'w') as f:
    #     f.write(str(counter)+'\n')
    #     for match in bestMatches:
    #         f.write(str(match[0])+' '+str(match[1])+'\n')
    #
    # with open("all_pairs_jaccard_standard_k_"+str(k)+".txt", 'w') as f:
    #     f.write(str(counter) + ' ' + str(numReads) + ' ' +
    #             str(time.clock()-timer) + '\n')
    #     for match in bestMatches:
    #         # f.write(str(match[0])+' '+str(match[1])+'\n')
    #         f.write(str(match)+'\n')

    for i in xrange(10, len(bestMatches)):
            with open("all_pairs_naive_jaccard_" + sim_measure + "_seqLength" +
                      str(i+1) + ".txt", 'w') as f:
                f.write(str(counter) + ' ' + str(numReads) + ' ' +
                        str(time.clock()-timer) + '\n')
                for match in bestMatches[i]:
                    f.write(str(match[0])+' '+str(match[1])+'\n')


# *************************************************************************** #
#                                                                             #
#                                     Main                                    #
#                                                                             #
# *************************************************************************** #
def main():
    """
    Main method of the program
    """
    totim = time.clock()

    # Parse command line options
    fasta_file, k, threshold, bands, rows, sim_measure, minhash_alg, \
        seed, log_file = optionParse()

    # n (number of hash functions = length of minhash signatures) is implicitly
    # given
    n = bands * rows

    if n % bands != 0:
        print "ERROR: The number of bands and rows do not go into n"
        sys.exit()

    with open(log_file, 'w') as log:

        #read_data_seq(fasta_file)
        LSH_run(fasta_file, bands, rows, n, k, seed, minhash_alg, log)
        sys.exit()

        # Compute all shingles in the data set
        shingles = computeAllShingles(fasta_file, k, log)

        # Read all reads from fasta file
        documents = readData(fasta_file, log)
        # for doc in readDataChunks(fasta_file, log):
        #     print doc.id

        # Compute all shingles in the data set
        # shingles = shingling(documents, k, log)

        # compareAllPairs(documents, k, sim_measure)
        # sys.exit()

        logprint(log, False, "Seed:", seed)
        if minhash_alg == 1:
            minhashingNew(documents, shingles, n, k, seed, log)
        elif minhash_alg == 2:
            minhashingOld(documents, shingles, n, k, seed, log)
        elif minhash_alg == 3:
            minhashingOld2(documents, shingles, n, k, seed, log)
        #minhashingOld3(documents, shingles, n, k, seed, log)
        LSH_old(documents, bands, rows, log)
        # LSH(documents, bands, rows, shingles, k, seed, log)

        # band_buckets = LSH_Old(documents, bands, rows)
        # LSH(documents, bands, rows, log)
        # LSH_Old(documents, bands, rows)

        logprint(log, True, "Total time used for LSH:",
                 (time.clock() - totim) / 60, "minutes")
        sys.exit()



        output_path = ("output_" + sim_measure + "/output_k_" + str(k) + "_b_"
                       + str(bands) + "_r_" + str(rows) + "/")
        output_file = ("output_k_" + str(k) + "_b_" + str(bands) + "_r_" +
                       str(rows) + ".txt")
        output = output_path + output_file
        # findSimilarPairsOld(band_buckets, threshold, k, totim, output,
        #                 len(documents)/2, sim_measure)
        findSimilarPairs(documents, threshold, k, output, len(documents)/2,
                         sim_measure)

        # bucketSize(band_buckets)

        print "Total time used:", time.clock() - totim / 60, "minutes"


if __name__ == '__main__':

    main()

    sys.exit(0)
