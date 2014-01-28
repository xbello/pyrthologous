"""Run all the functions with the genomes."""

import os

import blast
from .config import *
from .prepare import translate_fasta


if __name__ == "__main__":
    for pair in COMPARE:
        # Test if the files exists
        for specie in pair:
            if not os.path.isfile(os.path.join(BASE_PATH, GENOMES, specie)):
                raise IOError("File {0} doesn't exist".format(specie))

            # Translate both genomes only if they doesn't exist
            this_genome = os.path.join(BASE_PATH, GENOMES, specie)
            if not os.path.isfile(this_genome):
                translate_fasta(os.path.join(BASE_PATH, GENOMES, specie),
                                output_path=os.dirname(this_genome))

        # TODO: Reciprocal blast'em
        abs_paths = [os.path.join(BASE_PATH, OUTPUT, x) for x in pair]

        print [x for x in blast.reciprocal_blastp(abs_paths)]
        break


        # TODO: Select best matches and write'm out
        # TODO: Align the matches
        # TODO: Un-translate the sequences
        # TODO: Pack the results and clean the house
