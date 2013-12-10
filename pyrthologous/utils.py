"""Miscelaneous functions for the module."""

import itertools
import os
from Bio import SeqIO


def make_reciprocals(genomes):
    """Return a list of the permutations to generate from genomes."""

    return [x for x in itertools.combinations(genomes, 2)]


def normalize(fasta_file):
    """Yield a (id, head) for each line in a fasta file."""

    seq_number = 1
    if os.path.isfile(fasta_file):
        for rec in SeqIO.parse(fasta_file, "fasta"):
            yield seq_number, rec.name
            seq_number += 1


def reciprocal_genomes(src_dir, extension):
    """Return a list of the mutual searches to coduct."""

    genomes = select_genomes(src_dir, extension)

    return make_reciprocals(genomes)


def select_genomes(src_dir, extension):
    """Return the files in src_dir ending with extension."""

    if os.path.isdir(src_dir):
        return [x for x in os.listdir(src_dir) if x.endswith(extension)]


def tag_fasta(fasta, tag="NORM"):
    """Return a new name for a given fasta."""

    print fasta
    if os.path.isfile(fasta):
        new_fasta = list(os.path.split(os.path.abspath(fasta)))
        # Tag the new filename
        new_fasta[1] = new_fasta[1].split(".")[0] + "_" + tag + ".fas"

        return os.path.join(*new_fasta)
