#!/usr/bin/env python
"""Change the head/name of the fasta file to enter the BLAST commands.

- Names are unique by appending an int key.
- Save the old and new names in a file of equivalences.

"""

import os
from Bio import SeqIO, Seq

from .config import SEP
from .utils import tag_fasta


def process_file(fasta="", equivalences="", translate=False):
    """Return the status of the result of file processing."""

    if not os.path.isfile(fasta):
        return False

    new_fasta_handler = open(tag_fasta(fasta), "w")

    #new_fasta = list(os.path.split(os.path.abspath(fasta)))
    #new_fasta[1] = new_fasta[1].split(".")[0] + "_NORM.fas"
    #new_fasta_fn = os.path.join(new_fasta[0], new_fasta[1])
    #new_fasta_fp = open(new_fasta_fn, "w")

    equiv_file = open(equivalences, "w")

    n_seq = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        new_fasta_handler.write(">{0}\n".format(n_seq))
        equiv_file.write("{0}{1}{2}\n".format(rec.name, SEP, n_seq))
        n_seq += 1
        if translate:
            # Inter specie searches are better performed if a previous
            # translation of the sequence is being done:
            # http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&
            # PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#blastn

            new_fasta_handler.write("{0}\n".format(Seq.translate(rec.seq)))
        else:
            new_fasta_handler.write("{0}\n".format(rec.seq))

    equiv_file.close()
    new_fasta_handler.close()
    return True

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "fasta",
        help="Name of the file to normalize")
    parser.add_argument(
        "equivalences",
        help="Name of the file where equivalences will be writen")
    parser.add_argument(
        "-t", dest="translate", action="store_true",
        help="Spit out translated sequences instead of direct copies")

    args = parser.parse_args()

    print process_file(fasta=args.fasta,
                       equivalences=args.equivalences,
                       translate=args.translate)
