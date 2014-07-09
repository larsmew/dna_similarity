#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from __future__ import print_function

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "03/07/2014"
__version__ = "$Revision: 1.0"

from optparse import OptionParser
# import itertools
import sys
import time
import random
# import os


class Document(object):
    dna = ""
    id = 0
    isLeft = False  # If document comes from the left or right part of a read
    shingles = set()
    vector = []
    signature = []

    # Initializer
    def __init__(self, dna, id, isLeft):
        self.dna = dna
        self.id = id
        self.isLeft = isLeft


def createDocument(dna, id, isLeft):
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
                      default=5,
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
                      default=1/2,
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
    id = 1
    documents = []
    for read in reads:
        leftpart, rightpart = read[:len(read)/2], read[len(read)/2:]
        documents.append(createDocument(leftpart, id, True))
        documents.append(createDocument(rightpart, id, False))
        id += 1
    return documents


def shingling(docset, k):
    shingles = set()
    for doc in docset:
        docShingles = set()
        for i in range(len(doc.dna)-k):
            shingle = doc.dna[i:i+k]
            if shingle not in shingles:
                shingles.add(shingle)
            if shingle not in docShingles:
                docShingles.add(shingle)
        doc.shingles = docShingles
    return sorted(shingles)


def generateVectors(docset, shingles):
    for doc in docset:
        vector = []
        for shingle in shingles:
            if shingle in doc.shingles:
                vector.append(1)
            else:
                vector.append(0)
        doc.vector = vector
        print vector


def minhashing(docset, shingles, n):
    hashfuncs = []
    for i in range(n):
        h = range(len(shingles))
        random.shuffle(h)
        hashfuncs.append(h)
        print h,"\n"
    for doc in docset:
        #print docset[0].shingles
        signature = [None for i in range(n)]
        #print signature
        #print len(signature)
        for r in range(len(shingles)):
            if shingles[r] in doc.shingles:
                for i in range(len(hashfuncs)):
                    if signature[i] == None or signature[i] > hashfuncs[i][r]:
                        signature[i] = hashfuncs[i][r]
        #print signature


def main():
    """
    Main method of the program
    """
    totim = time.clock()

    # Parse command line options #
    fasta_file, k, n, threshold, bands, rows = optionParse()

    # Read all reads from fasta file #
    reads = readData(fasta_file)

    documents = createDocuments(reads)
    # for left,right in itertools.izip(leftDocs, rightDocs):
    #     print left.dna, left.id, left.isLeft
    #     print right.dna, right.id, right.isLeft

    shingles = shingling(documents, k)
    
    print shingles
    print len(shingles)

    minhashing(documents, shingles, n)

    print "Total time used:", time.clock() - totim


if __name__ == '__main__':

    main()

    sys.exit(0)
