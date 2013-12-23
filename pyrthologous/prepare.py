"""Prepare the files in the input to start the analysis."""

from .config import BASE_PATH, GENOMES
import os


def get_fasta_names(path=""):
    """Return all the filenames that can be compared."""
    if not path:
        path = os.path.join(BASE_PATH, GENOMES)

    return [z for x, y, z in os.walk(path)][0]
