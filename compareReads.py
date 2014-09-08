#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from __future__ import print_function

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/07/2014"
__version__ = "$Revision: 1.0"

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


class Document(object):
    """
    Defines the data structure of a document. Like dna sequence, id, left or
    right part of string, shingles related to the string, the minhash signature
    """
    dna = ""
    id = 0
    isLeft = False  # If document comes from the left or right part of a read
    shingles = set()
    signature = []
    similarTo = set()

    # Initializer
    def __init__(self, dna, id, isLeft):
        self.dna = dna
        self.id = id
        self.isLeft = isLeft


def createDocument(dna, id, isLeft):
    """
    Function to create document objects
    """
    document = Document(dna, id, isLeft)
    return document


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

    parser.add_option("-n", "--n_length",
                      metavar="<VALUE>",
                      type=int,
                      action="store",
                      dest="n",
                      help="set <VALUE> as length of minhash signatures.")

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
                      
    parser.add_option("-s", "--similarity_measure",
                      metavar="<VALUE>",
                      type=str,
                      default="naive",
                      action="store",
                      dest="s",
                      help="set <VALUE> as similairy measure to use.")

    (options, args) = parser.parse_args()

    return options.fasta_file, options.k, options.n, options.threshold,\
        options.bands, options.rows, options.s


def readData(fasta_file):
    """
    Extract the reads (DNA sequences) from the given fasta file
    """
    if fasta_file:
        reads = []
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
                    read += line.strip().upper()
            # Append the last read in the file to the list
            if read != "":
                reads.append(read)
        return reads
    else:
        print("ERROR: NO FASTA FILE GIVEN")


def createDocuments(reads):
    """
    Splits each read into two parts - left and right halfs - and creates
    the document objects.
    """
    id = 1  # Uses international id system for bookkeeping
    documents = []
    for read in reads:
        # Splits the string into two parts
        leftpart, rightpart = read[:len(read)/2], read[len(read)/2:]
        # Creates a document object for the left part of the read
        documents.append(createDocument(leftpart, id, True))
        # Creates a document object for the right part of the read
        documents.append(createDocument(rightpart, id, False))
        id += 1
    return documents


def shingling(docs, k):
    """
    Create k-shingles (substrings of length k) from each document
    """
    shingles = set()  # Contains all k-shingles in the dataset
    for doc in docs:
        docShingles = set()  # Contains the k-shingles in the given document
        for i in xrange(len(doc.dna)-k):
            # create k-shingle (substring) from the document
            shingle = doc.dna[i:i+k]
            # Add it to the set of all shingles
            shingles.add(shingle)
            # Add it to the set of shingles of the current document
            docShingles.add(shingle)
        # Store the document's shingles in the document object
        doc.shingles = docShingles
        # print docShingles,"\n"
    return shingles


def minhashing(docs, shingles, n):
    """
    Create minhash signatures using the shingles
    """
    # Create n different permutations (hash functions) of the shingles
    hashfuncs = []
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


def LSH(documents, bands, rows):
    """
    Use Locality-Sensitive Hashing (LSH) with the banding technique to
    hash similar elements to the same buckets.
    """
    band_buckets = []
    # buckets = dict([])
    index = 0
    counter = 0
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
            counter += 1
        print counter
        index += rows
        band_buckets.append(buckets)

    # for bucket in buckets:
    #     print "Bucket:", bucket
    #     for document in buckets[bucket]:
    #         print document.id,
    #     print
    #     for document in buckets[bucket]:
    #         print document.signature,
    #     print

    return band_buckets



#*****************************************************************************#
#                                                                             #
#                            Similarity checkers                              #
#                                                                             #
#*****************************************************************************#
def findSimilarPairs(band_buckets, t, totim, output, numReads, sim_measure):
    """
    Find candidate pairs that has a similarity above the threshold t
    """
    counter = 0
    counter2 = 0
    timer = time.clock()
    bestMatches = []
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
                                    if sim_measure == "naive":
                                        bestMatch = globalAlignment(doc1, doc2)
                                    elif sim_measure == "jaccard":
                                        bestMatch = jaccardSim(doc1, doc2)
                                else:
                                    if sim_measure == "naive":
                                        bestMatch = globalAlignment(doc2, doc1)
                                    elif sim_measure == "jaccard":
                                        bestMatch = jaccardSim(doc2, doc1)

                                bestMatches.append(bestMatch)

                                if counter2 % 500000 == 0:
                                    print "Processed", counter2, "pairs in", \
                                        "time:", time.clock() - timer
                            else:
                                print "Total time used (in secs):", \
                                    time.clock() - totim
                                bestMatches.sort(key=itemgetter(0),
                                                 reverse=False)
                                with open("output.txt", 'w') as f:
                                    counter = 0
                                    for match in bestMatches:
                                        f.write(str(counter) + "," + 
                                                str(match[0])+"\n")
                                        counter += 1
                                sys.exit()

    processing_time = time.clock() - timer
    print "Processed", counter2, "pairs in time:", processing_time
    print "{:,}".format(counter)
    print "{:,}".format(counter2)

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


def globalAlignment(doc1, doc2):
    # dna1, dna2 = doc1.dna, doc2.dna
    # print doc1.dna
    # print len(doc1.dna)
    # print doc2.dna
    # print len(doc2.dna)
    # print

    start = 0
    bestScore = (0,0,0)
    seqLength = len(doc1.dna)-start
    while seqLength > 10 and bestScore[0] != 1.0:
        #print seqLength, bestScore[1]
        matches = 0
        for i in xrange(seqLength):
            # print len(doc1.dna)-start
            if doc1.dna[i] == doc2.dna[i+start]:
                matches += 1
        #print bestScore
        score = matches / float(seqLength)
        if score > bestScore[0]:
            #print score, bestScore[0]
            #print seqLength, matches, bestScore[1]
            bestScore = (score, matches, seqLength)
        # if score > bestScore:
        #     bestScore = score
        start += 1
        seqLength = len(doc1.dna)-start
    return bestScore


def jaccardSim(doc1, doc2):
    intersection = doc1.shingles.intersection(doc2.shingles)
    if len(intersection) == 0:
        return 0
    union = doc1.shingles.union(doc2.shingles)
    return float(len(intersection)) / len(union)
    
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

    # For testing
    # print doc1.id, doc2.id
    # print doc1.dna
    # print doc2.dna
    # print counterA
    # print counterB
    # print intersection
    # print union
    # print float(intersection) / union

    # Check if jaccard similarity is above t
    # if float(intersection) / union >= t:
    #     doc1.similarTo.add(doc2)
    #     doc2.similarTo.add(doc1)


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


# test_longest_common_substring
def testLCS(doc1, doc2):
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


# NOT relevant method, just for playing around #
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


def main():
    """
    Main method of the program
    """
    totim = time.clock()

    # Parse command line options
    fasta_file, k, n, threshold, bands, rows, sim_measure = optionParse()

    if not n:
        n = bands*rows

    if bands*rows != n:
        print "ERROR: bands * rows do not equal n (number of hash functions)"
        sys.exit(0)

    # Read all reads from fasta file
    reads = readData(fasta_file)            

    documents = createDocuments(reads)
    # for doc in documents:
    #     print doc.dna, doc.id, doc.isLeft
    
    shingles = list(shingling(documents, k))
    # print shingles
    # print len(shingles)

    numReads = "{:,}".format(len(reads)).replace(",",".")
    print "Number of reads:", numReads
    print "Number of documents:", "{:,}".format(len(documents)).replace(",",".")
    numPairs = len(reads) * (len(reads)-1)
    strNumPairs = format(numPairs, ',d').replace(",",".")
    print "Number of pairs to process:", strNumPairs
    bestMatches = []
    counter = 0
    timer = time.clock()
    for i in xrange(0, len(documents)):
        for j in xrange(i+1, len(documents)):
            doc1, doc2 = documents[i], documents[j]
            # Check if doc1 and doc2 are candidate pairs
            if doc1.isLeft != doc2.isLeft and doc1.id != doc2.id:
                counter += 1
                if doc1.isLeft:
                    #bestMatch_naive = globalAlignment(doc1, doc2)
                    bestMatch_jaccard = jaccardSim(doc1, doc2)
                else:
                    #bestMatch_naive = globalAlignment(doc2, doc1)
                    bestMatch_jaccard = jaccardSim(doc2, doc1)
                #bestMatches.append((bestMatch_naive[0], bestMatch_jaccard))
                bestMatches.append(bestMatch_jaccard)
                
                if counter % 200000 == 0:
                    print "Processed", format(counter, ',d').replace(",","."), \
                          "of", strNumPairs, "pairs in time (minutes):", \
                          (time.clock() - timer) / 60
                    
    # bestMatches.sort(key=itemgetter(0), reverse=False)
    # with open("all_pairs_naive.txt", 'w') as f:
    #     f.write(str(counter)+'\n')
    #     for match in bestMatches:
    #         f.write(str(match[0])+' '+str(match[1])+'\n')

    #bestMatches.sort(key=itemgetter(1), reverse=False)
    with open("all_pairs_jaccard_normal.txt", 'w') as f:
        f.write(str(counter)+' '+str(numReads)+' '+str(time.clock()-timer)+'\n')
        for match in bestMatches:
            #f.write(str(match[0])+' '+str(match[1])+'\n')
            f.write(str(match)+'\n')

    # minhashing(documents, shingles, n)
    #
    # band_buckets = LSH(documents, bands, rows)
    #
    # print "Time used:", time.clock() - totim
    #
    # output_path = ("output_"+sim_measure+"/output_k_" + str(k) + "_b_" +
    #                str(bands) + "_r_" + str(rows) + "/")
    # output_file = "output_k_"+str(k)+"_b_"+str(bands)+"_r_"+str(rows)+".txt"
    # output = output_path + output_file
    # numReads = len(reads)
    # findSimilarPairs(band_buckets, threshold, totim, output, numReads,
    #                  sim_measure)

    # bucketSize(band_buckets)

    print "Total time used (in secs):", time.clock() - totim


if __name__ == '__main__':

    main()

    sys.exit(0)
