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


def set_pairs(pair_line):
    """Return a tuple with the filenames of the files to be compared.
    
    Files to be compared are strings with the names of the fasta files joined
    with '_', valid formats are::
	    
      1_vs_2
      1_to_2
      1_2
      
    When the target files are 1.fas and 2.fas."""

    pairs = [x + ".fas" for x in pair_line.split("_")]

    return (pairs[0], pairs[-1])

def translate_fasta(fasta_path, output_path):
    """Return True if can make a translate file from a file fasta."""

    records = SeqIO.parse(fasta_path, "fasta")
    # Create a new fasta file to save the output
    output = open(
        os.path.join(output_path, os.path.basename(fasta_path)), "w")

    SeqIO.write(translate(records), output, "fasta")

    output.close()

    return True
