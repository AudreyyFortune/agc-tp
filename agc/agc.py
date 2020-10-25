#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Fortune Audrey"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Fortune Audrey"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Fortune Audrey"
__email__ = "audefortune@free.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
	"""
	Read a fasta file.
	Parameters: - amplicon_file (.fasta)
				- minseqlen (minimum length of sequences)
	Return: a generator of sequences of size
	"""
	with gzip.open(amplicon_file, "rb") as filin:
        for seq in filin:
            if len(seq) >= minseqlen:
                yield seq.strip()
            continue


def dereplication_fullength(amplicon_file, minseqlen, mincount):
	"""
	Determine unique sequences that occurs.
	Parameters: - amplicon_file (.fasta)
				- minseqlen (minimum length of sequences)
				- mincount (minimum counting)
	Return: a single sequence generator with occurrence [sequence, count]
	"""
	dict_occ = {}
	for seq in read_fasta(amplicon_file, minseqlen):
        if seq not in dict_occ.keys():    
            dict_occ[seq] = 0
        dict_occ[seq] += 1
        #dict_occ[seq] += 1
    new_dict_occ = dict(sorted(dict_occ.items(), key = lambda item : item[1], reverse = True))
    for seq, count in new_dict_occ.items():
        if count >= mincount:
            yield [seq, count]
	

def get_chunk(seq, chunk_size):
	"""
	To get a sub-sequences of size not overlapping.
	Parameters: - seq (sequence in fastq)
				- chunk_size (segment length)
	Return: a chunk 
	"""
	chunk_list = []
	for i in range(0, len(sequence), chunk_size):
		if chunk_list:
			if len(sequence[i:i + chunk_size]) != len(chunk_list[0]):
				break
		chunk_list.append(sequence[i:i + chunk_size])
	if len(chunk_list) > 5:  # >=4
		return chunk_list
		

def cut_kmer(seq, kmer_size):
	"""
	Cuts the sequence in kmer with a kmer size.
	Parameters: - seq (sequence in fastq)
				- kmer_size (size of kmer)
	Return: a generator of sequence's kmers.
	"""
	for i in range(len(seq) - kmer_size + 1):
		yield seq[i: i + kmer_size]


def get_unique_kmer(kmer_dict, seq, id_seq, kmer_size):
	"""
	Builds a k-mer dictionary for a giver sequence.
	Parameters: - kmer_dict (dictionnary of k-mer)
				- seq (sequence in fastq)
				- id_seq (id of the sequence)
				- kmer_size (size of kmer)
	Return: 
	"""
	for sequence in cut_kmer(seq, kmer_size):
		if sequence not in kmer_dict:
			kmer_dict[sequence] = [id_seq]
		elif id_seq not in kmer_dict[sequence]:
			kmer_dict[sequence].append(id_seq)
	return kmer_dict
	

def search_mates(kmer_dict, seq, kmer_size):
	"""
	"""
	pass
	

def get_identity(alignment_list):
	"""
	Calculates the percentage of identity between the two sequences.
	Parameters: - alignment_list (list of alignment)
	Return: the identity between the two sequences.
	
	"""
	compte = 0
	for i in range(len(alignment_list[0])):
		if alignment_list[0][i] == alignment_list[1][i]:
			compte += 1
	identity = (compte / len(alignment_list[0])) * 100
	return identity


def detect_chimera(perc_identity_matrix):
	"""
	"""
	pass


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""
	"""
	pass
	

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""
	"""
	pass


def fill(text, width = 80):
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
	"""
	"""
	pass



#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()