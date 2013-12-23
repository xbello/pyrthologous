"""Prepare the files in the input to start the analysis."""

from Bio import SeqIO

from .config import BASE_PATH, GENOMES, OUTPUT
import os


def get_fasta_names(path=""):
    """Return all the filenames that can be compared."""
    if not path:
        path = os.path.join(BASE_PATH, GENOMES)

    return [z for x, y, z in os.walk(path)][0]


def translate_fasta(fasta_path, output_path):
    """Make a translate file from a file fasta."""

    records = SeqIO.parse(fasta_path, "fasta")


