#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from __future__ import print_function

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/07/2014"
__version__ = "$Revision: 1.0"

from optparse import OptionParser
from collections import Counter
# import itertools
import sys
import time
import random
# import os


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
                      default=100,
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

    (options, args) = parser.parse_args()

    return options.fasta_file, options.k, options.n, options.threshold,\
        options.bands, options.rows


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
        for i in range(len(doc.dna)-k):
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
    for i in range(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        # print h,"\n"

    # Create minhash signatures as described in chapter 3 of the book Massive
    # Data Mining
    # Find signature for each document
    for doc in docs:
        signature = [None for i in range(n)]
        # For each row in the 'character matrix'
        for r in range(len(shingles)):
            # If the shingle is in the document, then
            if shingles[r] in doc.shingles:
                # Find the 'first' shingle relative to each permutation
                for i in range(len(hashfuncs)):
                    if signature[i] is None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        doc.signature = signature
        # print signature


def LSH(documents, bands, rows):
    band_buckets = []
    # buckets = dict([])
    index = 0
    counter = 0
    for b in range(bands):
        buckets = dict([])
        for i in range(len(documents)):
            col = documents[i].signature[index:index+rows]
            key = ''.join(map(str, col))
            if key in buckets:
                buckets[key].append(documents[i])
            else:
                buckets[key] = [documents[i]]
            counter += 1
        print counter
        index += rows
        band_buckets.append(buckets)

    # for bucket in buckets:
    #     print "Bucket:", bucket
    #     for document in buckets[bucket]:
    #         print document.id,
    #     print

    return band_buckets


def findSimilarPairs(band_buckets, t):
    counter = 0
    counter2 = 0
    for buckets in band_buckets:
        for bucket in buckets:
            if len(buckets[bucket]) > 1:
                for i in range(len(buckets[bucket])):
                    for j in range(i+1, len(buckets[bucket])):
                        doc1, doc2 = buckets[bucket][i], buckets[bucket][j]
                        counter += 1
                        if doc1.id != doc2.id and doc1.isLeft != doc2.isLeft:
                            counterA = Counter(doc1.signature)
                            counterB = Counter(doc2.signature)
                            intersection = sum((counterA & counterB).values())
                            union = sum((counterA | counterB).values())
                            if counter2 == 0:
                                print counterA
                                print counterB
                                print intersection
                                print union
                                print float(intersection) / union
                            if float(intersection) / union >= t:
                                # print doc1.id, doc2.id
                                # print doc1.isLeft, doc2.isLeft
                                # print intersection, union
                                # print float(intersection) / union
                                # print
                                counter2 += 1
                                doc1.similarTo.add(doc2)
                                doc2.similarTo.add(doc1)
    print "{:,}".format(counter)
    print "{:,}".format(counter2)


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

    # Parse command line options #
    fasta_file, k, n, threshold, bands, rows = optionParse()

    if bands*rows != n:
        print "ERROR: bands * rows not equal to n (number of hash functions)"

    # Read all reads from fasta file #
    reads = readData(fasta_file)

    documents = createDocuments(reads)
    # for doc in documents:
    #     print doc.dna, doc.id, doc.isLeft

    shingles = list(shingling(documents, k))

    # print shingles
    # print len(shingles)

    minhashing(documents, shingles, n)

    band_buckets = LSH(documents, bands, rows)

    print "Time used:", time.clock() - totim

    findSimilarPairs(band_buckets, threshold)
    # bucketSize(band_buckets)

    print "Total time used (in secs):", time.clock() - totim


if __name__ == '__main__':

    main()

    sys.exit(0)
