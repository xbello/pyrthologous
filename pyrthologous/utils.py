"""Miscelaneous functions for the module."""

import itertools
import os
from .config import BASE_PATH


def make_reciprocals(genomes):
    """Return a list of the permutations to generate from genomes."""

    return [x for x in itertools.combinations(genomes, 2)]


def reciprocal_genomes(src_dir, extension):
    """Return a list of the mutual searches to coduct."""

    genomes = select_genomes(src_dir, extension)

    return make_reciprocals(genomes)


def select_genomes(src_dir, extension):
    """Return the files in src_dir relative to BASE_PATH with extension."""

    genome_src = os.path.join(BASE_PATH, src_dir)

    if os.path.isdir(genome_src):
        return [x for x in os.listdir(genome_src) if x.endswith(extension)]


def tag_fasta(fasta, tag="NORM"):
    """Return a new name for a given fasta."""

    if os.path.isfile(fasta):
        new_fasta = list(os.path.split(os.path.abspath(fasta)))
        # Tag the new filename
        new_fasta[1] = new_fasta[1].split(".")[0] + "_" + tag + ".fas"

        return os.path.join(*new_fasta)
