#!/usr/bin/env python
import os
from Bio import SeqIO

def process_file(fasta="", equivalences=""):
    """filename, filename -> bool
    """
    if not os.path.isfile(fasta):
        return False

    new_fasta = list(os.path.split(os.path.abspath(fasta)))
    new_fasta[1] = new_fasta[1].split(".")[0] + "_NORM.fas"
    new_fasta_fn = os.path.join(new_fasta[0], new_fasta[1])
    new_fasta_fp = open(new_fasta_fn, "w")

    equiv_file = open(equivalences, "w") 

    n_seq = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        new_fasta_fp.write(">{0}\n".format(n_seq))
        equiv_file.write("{0},{1}\n".format(rec.name, n_seq))
        n_seq += 1
        new_fasta_fp.write("{0}\n".format(rec.seq))

    equiv_file.close()
    new_fasta_fp.close()
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
