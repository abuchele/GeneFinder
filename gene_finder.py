# -*- coding: utf-8 -*-
"""
GENE FINDER

@author: Anna Buchele

Takes dna strand and finds all open reading frames in both directions,
    compares their length to the longest open reading frame found in 
    a randomly shuffled version of the same dna strand, and translates
    into the corresponding amino acids those that are longer than the
    threshold length of the longest open reading frame found in the 
    shuffled dna. 
Returns the translated list of open reading frames.

"""

import random
from amino_acids import aa, codons, aa_table 
from load import load_seq
import numpy as np
from operator import add


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a 			string
        returns: the complementary nucleotide"""
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'
    else:
        return False


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the 		specfied DNA sequence
		dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as 			a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    reverse=dna[::-1]
    revdna=''
    for x in reverse:
    	revnuc=get_complement(x)
    	revdna+=revnuc
    return revdna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stopcodons= ['TAA','TAG','TGA']
    i = 0
    while i < len(dna):
    	if dna[i:i+3] in stopcodons:
    		RORF=dna[0:i]
    		return RORF
    	else:
    		i=i+3
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    remainingdna=dna
    i=0
    oneframeORFs=[]
    while i<len(dna):
        if dna[i:i+3]=='ATG':
            remainingdna=dna[i:]
            orfs= rest_of_ORF(remainingdna)
            oneframeORFs.append(orfs)
            i=i+len(orfs)
        else:
            i=i+3
    return oneframeORFs

        
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORFS=[]
    for x in range (0,3):
        dnaframe=dna[x:]
        ORFSnew=find_all_ORFs_oneframe(dnaframe)
        if len(ORFSnew)==0:
        	continue
        ORFS+=ORFSnew
    return ORFS
        

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    revcompdna=get_reverse_complement(dna)
    ORFS1= find_all_ORFs(dna)
    ORFS2= find_all_ORFs(revcompdna)
    ORFS1+=ORFS2
    return ORFS1

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFS=find_all_ORFs_both_strands(dna)
    toplength=0
    n=len(ORFS)
    for x in range (0,n):
        thisone=ORFS[x]
        thislength=len(thisone)
        if thislength > toplength:
            toplength= thislength
            longestORF= ORFS[x]
    return longestORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    for x in range (0,num_trials):
        shuffle= shuffle_string(dna)
    maxlengthORF= longest_ORF(shuffle)
    return maxlengthORF

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    aminos=''
    length=len(dna)
    for x in range (0,length,3):
        code=dna[x:x+3]
        if code=='ATG':
            aminos= aminos+'M'
        elif code== 'TTT':
            aminos+='F'
        elif code== 'TTC':
            aminos+='F'
        elif code== 'TTA':
            aminos+='L'
        elif code== 'TTG':
            aminos+='L'
        elif code== 'CTT':
            aminos+='L'
        elif code== 'CTG':
            aminos+='L'
        elif code== 'CTC':
            aminos+='L'
        elif code== 'CTA':
            aminos+='L'
        elif code== 'ATA':
            aminos+='I'
        elif code== 'ATT':
            aminos+='I'
        elif code== 'ATC':
            aminos+='I'
        elif code== 'GTT':
            aminos+='V'
        elif code== 'GTC':
            aminos+='V'
        elif code== 'GTA':
            aminos+='V'
        elif code== 'GTG':
            aminos+='V'
        elif code== 'TCT':
            aminos+='S'
        elif code== 'TCC':
            aminos+='S'
        elif code== 'TCA':
            aminos+='S'
        elif code== 'TCG':
            aminos+='S'
        elif code== 'CCT':
            aminos+='P'
        elif code== 'CCC':
            aminos+='P'
        elif code== 'CCA':
            aminos+='P'
        elif code== 'CCG':
            aminos+='P'
        elif code== 'GCT':
            aminos+='A'
        elif code== 'GCC':
            aminos+='A'
        elif code== 'GCA':
            aminos+='A'
        elif code== 'GCG':
            aminos+='A'
        elif code== 'TAT':
            aminos+='Y'
        elif code== 'TAC':
            aminos+='Y'
        elif code== 'CAT':
            aminos+='H'
        elif code== 'CAC':
            aminos+='H'
        elif code== 'CAA':
            aminos+='Q'
        elif code== 'CAG':
            aminos+='Q'
        elif code== 'CGT':
            aminos+='R'
        elif code== 'CGC':
            aminos+='R'
        elif code== 'CGA':
            aminos+='R'
        elif code== 'CGG':
            aminos+='R'
        elif code== 'GGT':
            aminos+='G'
        elif code== 'GGC':
            aminos+='G'
        elif code== 'GGA':
            aminos+='G'
        elif code== 'GGG':
            aminos+='G'
        elif code== 'ATT':
            aminos+='N'
        elif code== 'AAC':
            aminos+='N'
        elif code== 'AAA':
            aminos+='K'
        elif code== 'AAG':
            aminos+='K'
        elif code== 'GAT':
            aminos+='D'
        elif code== 'GAC':
            aminos+='D'
        elif code== 'GAA':
            aminos+='E'
        elif code== 'GAG':
            aminos+='E'
        elif code== 'TGT':
            aminos+='C'
        elif code== 'TGC':
            aminos+='C'
        elif code== 'TGG':
            aminos+='W'
        elif code== 'AGT':
            aminos+='S'
        elif code== 'AGC':
            aminos+='S'
        elif code== 'AGA':
            aminos+='R'
        elif code== 'AGG':
            aminos+='R'
        elif code== 'TAA':
            aminos+='*'
        elif code== 'TAG':
            aminos+='*'
        elif code== 'TGA':
            aminos+='*'
    return aminos


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    above_thresh_ORFs=[]
    amino_seq=''
    threshold = longest_ORF_noncoding(dna, 1500)
    ALLORFS= find_all_ORFs_both_strands(dna)
    NumORFS=len(ALLORFS)
    for x in range (0, NumORFS):
        thisorf=ALLORFS[x]
        ORFlen= len(thisorf)
        if ORFlen > len(threshold):
            above_thresh_ORFs.append(thisorf)
    count=0
    NumATORFS=len(above_thresh_ORFs)
    for x in range (0,NumATORFS):
        aminos= coding_strand_to_AA(above_thresh_ORFs[x])
        if count==0:
            finaldna=list([aminos])
        else:
            finaldna+=[aminos]
        count=count+1
    return finaldna

from load import load_seq
dna = load_seq("./data/X73525.fa")

print gene_finder(dna)

# Doctest:
# if __name__ == "__main__":
#     import doctest
#     doctest.testmod()
