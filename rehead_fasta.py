#!/usr/bin/env python
import os
from Bio import SeqIO

def process_file(fasta="", equivalences=""):
    """filename, filename -> bool
    """
    if not os.path.isfile(fasta):
        return False

    equiv_file = open(equivalences, "w") 
    return True

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta",
        help="Name of the file to normalize")
    parser.add_argument("equivalences",
        help="Name of the file where equivalences will be writen")

    args = parser.parse_args()

    print process_file(fasta = args.fasta,
        equivalences = args.equivalences)
