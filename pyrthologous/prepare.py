"""Prepare the files in the input to start the analysis."""

from Bio import SeqIO

from .config import BASE_PATH, GENOMES
from .utils import translate
import os


def get_fasta_names(path=""):
    """Return all the filenames that can be compared."""
    if not path:
        path = os.path.join(BASE_PATH, GENOMES)

    return [x[2] for x in os.walk(path)][0]


def translate_fasta(fasta_path, output_path):
    """Return True if can make a translate file from a file fasta."""

    records = SeqIO.parse(fasta_path, "fasta")
    # Create a new fasta file to save the output
    output = open(
        os.path.join(output_path, os.path.basename(fasta_path)), "w")

    SeqIO.write(translate(records), output, "fasta")

    output.close()

    return True
