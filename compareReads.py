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
# import os
# import math
# import numpy
# import nwalign as nw


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
    shingles = set()
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
                      default="./Data/reads.fa",
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
                      default=4,
                      action="store",
                      dest="bands",
                      help="set <VALUE> as the number of bands for LSH.")

    parser.add_option("-r", "--rows",
                      metavar="<Value>",
                      type=int,
                      default=3,
                      action="store",
                      dest="rows",
                      help="set <VALUE> as the number of rows for LSH.")

    parser.add_option("-m", "--similarity_measure",
                      metavar="<VALUE>",
                      type=str,
                      default="naive",
                      action="store",
                      dest="m",
                      help="set <VALUE> as similairy measure to use.")

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
        options.bands, options.rows, options.m, options.s, options.log


def memory_usage_resource():
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def readDataOld(fasta_file):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    if fasta_file:
        print "Processing fasta file..."
        tim = time.clock()
        reads = []
        seqs = 0
        with open(fasta_file, "rU") as fasta_file:
            read = ""
            for line in fasta_file:
                # If line starts with ">", which indicates end of a sequence
                # Then, append it to list of reads
                if line.startswith(">"):
                    if read != "":
                        reads.append(read)
                        read = ""
                # Concatenate multi-line sequences into one string
                else:
                    seqs += 1
                    read += line.strip().upper()
            # Append the last read in the file to the list
            if read != "":
                reads.append(read)

        print "Read", format(seqs, ',d'), "sequences in time:", \
              (time.clock() - tim) / 60, "minutes"
        return reads
    else:
        print("ERROR: NO FASTA FILE GIVEN")


def readData(fasta_file, log):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    if fasta_file:
        print "Processing fasta file..."
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
                        # Creates a document object for the left part of the
                        # read
                        if id == 0:
                            print leftpart
                        documents.append(createDocument(leftpart, id, True))
                        id += 1
                        # Creates a document object for the right part of the
                        # read
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

        print "Read", format(seqs, ',d'), "sequences in", \
              (time.clock() - tim) / 60, "minutes"
        log.write("Read " + str(format(seqs, ',d')) + " sequences in " \
              + str((time.clock() - tim) / 60) + " minutes")
        print memory_usage_resource()
        return documents
    else:
        print("ERROR: NO FASTA FILE GIVEN")
        sys.exit()


# *************************************************************************** #
#                                                                             #
#                          Locality Sensitive Hashing                         #
#                                                                             #
# *************************************************************************** #
def createDocuments(reads):
    """
    Splits each read into two parts - left and right halfs - and creates
    the document objects.
    """
    tim = time.clock()
    print "Creating documents..."
    id = 1  # Uses id system for bookkeeping
    documents = []
    for read in reads:
        # Splits the string into two parts
        leftpart, rightpart = read[:len(read)/2], read[len(read)/2:]
        # Creates a document object for the left part of the read
        documents.append(createDocument(leftpart, id, True))
        # Creates a document object for the right part of the read
        documents.append(createDocument(rightpart, id, False))
        id += 1
    print "Created documents in", (time.clock() - tim) / 60, "minutes"
    return documents


# Obsolete - faster, but too high memory requirements
def shinglingOld(docs, k):
    """
    Create k-shingles (substrings of length k) from each document
    """
    tim = time.clock()
    print "Shingling..."
    shingles = set()  # Contains all k-shingles in the dataset
    count = 0
    for doc in docs:
        count += 1
        docShingles = set()  # Contains the k-shingles in the given document
        counterShingles = []
        for i in xrange(len(doc.dna)-k+1):
            # create k-shingle (substring) from the document
            shingle = doc.dna[i:i+k]
            # Add it to the set of all shingles
            shingles.add(shingle)
            # Add it to the set of shingles of the current document
            docShingles.add(shingle)
            # counterShingles.append(shingle)

        # doc.counterShingles = Counter(counterShingles)
        # Store the document's shingles in the document object
        doc.shingles = docShingles

        # print docShingles,"\n"
        if count % 500000 == 0:
            print "Processed", count, "of 27.057.600 documents"
    print "Finished shingling in", (time.clock() - tim) / 60, "minutes"
    return shingles


def getDocShingles(doc, k):
    shingles = set()
    for i in xrange(len(doc.dna)-k+1):
        # create k-shingle (substring) from the document
        shingle = doc.dna[i:i+k]
        # Add it to the set of all shingles
        shingles.add(shingle)
    return shingles


def shingling(docs, k, log):
    tim = time.clock()
    print "Shingling..."
    shingles = set()  # Contains all k-shingles in the dataset
    # shinglesDict = dict()
    for doc in docs:
        for i in xrange(len(doc.dna)-k+1):
            # create k-shingle (substring) from the document
            shingle = doc.dna[i:i+k]
            # Add it to the set of all shingles
            shingles.add(shingle)
            # if shingle in shinglesDict:
            #     shinglesDict[shingle].add(doc.id)
            # else:
            #     shinglesDict[shingle] = set([doc.id])
    print "Finished shingling in", (time.clock() - tim) / 60, "minutes"
    return list(shingles)


def minhashing(docs, shingles, n, k, seed, log):
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
    print "Minhashing..."
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    count = 0
    for doc in docs:
        count += 1
        docShingles = getDocShingles(doc, k)
        signature = [None for i in xrange(n)]
        # For each row in the 'character matrix'
        occupied = 0
        #while occupied < len(signature):
        for r in xrange(len(shingles)):
            if occupied == len(signature):
                break
            for i in xrange(n):
                if signature[i] is None:
                    if shingles[hashfuncs[i][r]] in docShingles:
                        signature[i] = r
                        occupied += 1
        doc.signature = signature
        # print signature
        if count % 100 == 0:
            print "Processed", count, "documents in", (time.clock() - tim) \
                  / 60, "minutes"

    print "Finished minhashing in", (time.clock() - tim) / 60, "minutes"


def minhashingOld(docs, shingles, n, k, seed, log):
    """
    Create minhash signatures using the shingles
    """
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    random.seed(seed)
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"
    print "Computed hashfunctions"

    tim = time.clock()
    print "Minhashing..."
    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    count = 0
    for doc in docs:
        count += 1
        docShingles = getDocShingles(doc, k)
        signature = [None for i in xrange(n)]
        # For each row in the 'character matrix'
        for r in xrange(len(shingles)):
            # If the shingle is in the document, then
            # if doc.id in shinglesDict[shingles[r]]:
            if shingles[r] in docShingles:
                # Find the 'first' shingle relative to each permutation
                for i in xrange(n):
                    if signature[i] is None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        # print sum(signature)
        doc.signature = signature
        # print signature
        if count % 100 == 0:
            print "Processed", count, "documents in", (time.clock() - tim) \
                  / 60, "minutes"


    print "Finished minhashing in", (time.clock() - tim) / 60, "minutes"


def minhashingOld_MemHungry(docs, shingles, n, seed):
    """
    Create minhash signatures using the shingles
    """
    tim = time.clock()
    print "Minhashing..."
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
    random.seed(seed)
    for i in xrange(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"

    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    for doc in docs:
        signature = [None for i in xrange(n)]
        # For each row in the 'character matrix'
        for r in xrange(len(shingles)):
            # If the shingle is in the document, then
            if shingles[r] in doc.shingles:
                # Find the 'first' shingle relative to each permutation
                for i in xrange(n):
                    if signature[i] is None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        doc.signature = signature
        # print signature
    print "Finished minhashing in", (time.clock() - tim) / 60, "minutes"


def LSH(documents, bands, rows, log):
    """
    Use Locality-Sensitive Hashing (LSH) with the banding technique to
    hash similar elements to the same buckets.
    """
    tim = time.clock()
    print "Running LSH..."
    index = 0
    for b in xrange(bands):
        # For each band, create a bucket array with signature as key
        buckets = dict([])
        for doc in documents:
            # Obtain sub-signature of length rows
            col = doc.signature[index:index+rows]
            # Convert to string
            key = ''.join(map(str, col))
            # Place the document in a bucket
            if key in buckets:
                buckets[key].append(doc)
            else:
                buckets[key] = [doc]
        index += rows

        for bucket in buckets:
            for i in xrange(len(buckets[bucket])):
                doc1 = buckets[bucket][i]
                for j in xrange(i+1, len(buckets[bucket])):
                    doc2 = buckets[bucket][j]
                    if doc1.isLeft and not doc2.isLeft:
                        if doc1.id + 1 != doc2.id:
                            doc1.similarTo.add(doc2)
                    if not doc1.isLeft and doc2.isLeft:
                        if doc1.id != doc2.id + 1:
                            doc1.similarTo.add(doc2)
                    # if doc1.isLeft != doc2.isLeft:
                    #     doc1.similarTo.add(doc2)

        print "Number of buckets in band", b+1, ":", len(buckets)
        numPairs = 0
        for bucket in buckets:
            inBucket = buckets[bucket]
            numPairs += len(inBucket) * (len(inBucket)-1) / 2
        print "Number of candidate pairs in band", b+1, ":", numPairs

    count = 0
    for doc in documents:
        count += len(doc.similarTo)
    print "Number of unique candidate pairs", count

    print "Finished LSH in", (time.clock() - tim) / 60, "minutes"
    return None

def LSH_Old(documents, bands, rows):
    """
    Use Locality-Sensitive Hashing (LSH) with the banding technique to
    hash similar elements to the same buckets.
    """
    tim = time.clock()
    print "Running LSH..."
    band_buckets = []
    index = 0
    for b in xrange(bands):
        # For each band, create a bucket array with signature as key
        buckets = dict([])
        for doc in documents:
            # Obtain sub-signature of length rows
            col = doc.signature[index:index+rows]
            # Convert to string
            key = ''.join(map(str, col))
            # Place the document in a bucket
            if key in buckets:
                buckets[key].append(doc)
            else:
                buckets[key] = [doc]
        index += rows
        band_buckets.append(buckets)

    for buckets in enumerate(band_buckets):
        print "Number of buckets in band", buckets[0]+1, ":", len(buckets[1])
        numPairs = 0
        for bucket in buckets[1]:
            inBucket = buckets[1][bucket]
            numPairs += len(inBucket) * (len(inBucket)-1) / 2
        print "Number of candidate pairs in band", buckets[0]+1, ":", numPairs

    # for buckets in band_buckets:
    #     for bucket in buckets:
    #         print "Bucket:", bucket
    #         for document in buckets[bucket]:
    #             print document.id,
    #         print
    #         for document in buckets[bucket]:
    #             print document.dna,
    #         print
    # sys.exit()

    for buckets in band_buckets:
        for bucket in buckets:
            for i in xrange(len(buckets[bucket])):
                pairs = set()
                doc1 = buckets[bucket][i]
                for j in xrange(i+1, len(buckets[bucket])):
                    doc2 = buckets[bucket][j]
                    if doc1.isLeft != doc2.isLeft: # and doc1.id != doc2.id:
                        pairs.add(doc2)
                doc1.similarTo = doc1.similarTo.union(pairs)

    count = 0
    for doc in documents:
        #print [i.id for i in doc.similarTo]
        count += len(doc.similarTo)
    print "Number of unique candidate pairs", count

    print "Finished LSH in", (time.clock() - tim) / 60, "minutes"
    return band_buckets


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
    fasta_file, k, threshold, bands, rows, sim_measure, seed, log_file \
        = optionParse()

    # n (number of hash functions = length of minhash signatures) is implicitly
    # given
    n = bands * rows

    if n % bands != 0:
        print "ERROR: The number of bands and rows do not go into n"
        sys.exit()

    with open(log_file, 'w') as log:
        # Read all reads from fasta file
        documents = readData(fasta_file, log)

        # documents = createDocuments(reads)
        # for doc in documents:
        #     print doc.dna, doc.id, doc.isLeft

        shingles = shingling(documents, k, log)
        print "Number of shingles:", len(shingles)
        # shingles = list(shinglingOld(documents, k))
        # print shingles
        # print len(shingles)

        # compareAllPairs(documents, k, sim_measure)
        # sys.exit()

        print "Seed:", seed
        minhashing(documents, shingles, n, k, seed, log)
        # minhashingOld(documents, shingles2, n, seed)

        # band_buckets = LSH_Old(documents, bands, rows)
        LSH(documents, bands, rows, log)
        # LSH_Old(documents, bands, rows)

        print "Total time used for LSH:", (time.clock() - totim) / 60, "min"
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

        print "Total time used (in secs):", time.clock() - totim


if __name__ == '__main__':

    main()

    sys.exit(0)
