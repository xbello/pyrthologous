"""Miscelaneous functions for the module."""

import itertools
import os
from .config import BASE_PATH


def make_reciprocals(genomes):
    """Return a list of the permutations to generate from genomes."""

    return [x for x in itertools.combinations(genomes, 2)]


def select_genomes(src_dir, extension):
    """Return the files in src_dir relative to BASE_PATH with extension."""

    genome_src = os.path.join(BASE_PATH, src_dir)

    if os.path.isdir(genome_src):
        return [x for x in os.listdir(genome_src) if x.endswith(extension)]
