"""Prepare the files in the input to start the analysis."""

import logging
import os
from itertools import chain
from Bio import SeqIO
from Bio.SeqIO import FastaIO

from .utils import translate

logger = logging.getLogger(__name__)


def check_compare(path="", compare=(), config=None):
    """Return True if all the genomes in compare exists in path."""
    if not path:
        path = os.path.join(config.BASE_PATH, config.GENOMES)

    if not compare:
        compare = config.COMPARE

    available_genomes = get_fasta_names(path=path)

    for first, second in compare:
        assert first in available_genomes
        assert second in available_genomes

    return True


def get_fasta_names(path="", config=None):
    """Return all the filenames that can be compared."""
    if not path:
        path = os.path.join(config.BASE_PATH, config.GENOMES)

    return [x[2] for x in os.walk(path)][0]


def translate_all_fastas(compare=(), output_path="", config=None):
    """Translate all fasta files in compare and output them in output_path."""
    if not compare:
        compare = config.COMPARE
    if not output_path:
        output_path = os.path.join(config.BASE_PATH, config.OUTPUT)

    # This get the filenames one and only one time
    filenames = set(chain.from_iterable(compare))

    logger.info("Creating new translated fasta in {0}".format(output_path))
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    for filename in filenames:
        translate_fasta(
            os.path.join(config.BASE_PATH, config.GENOMES, filename),
            output_path=output_path)

    return True


def translate_fasta(fasta_path, output_path="", config=None):
    """Return True if can make a translate file from a file fasta."""

    if not output_path:
        output_path = os.path.join(config.BASE_PATH, config.OUTPUT)

    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    records = FastaIO.FastaIterator(open(fasta_path))
    # Create a new fasta file to save the output
    with open(os.path.join(
            output_path, os.path.basename(fasta_path)), "w") as output:
        SeqIO.write(translate(records), output, "fasta")

    logger.info("Tranlated {0}".format(fasta_path))

    return True
