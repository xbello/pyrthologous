"""Run all the functions with the genomes."""

import os

from .config import *
from .prepare import translate_fasta


if __name__ == "__main__":
    for pair in COMPARE:
        # Test if the files exists
        for specie in pair:
            if not os.path.isfile(os.path.join(BASE_PATH, GENOMES, specie)):
                print os.path.join(BASE_PATH, GENOMES, specie)
                raise IOError("File {0} doesn't exist".format(specie))

            # Translate both genomes
            translate_fasta(os.path.join(BASE_PATH, GENOMES, specie),
                            output_path=os.path.join(BASE_PATH, OUTPUT))

        # TODO: Reciprocal blast'em
        # TODO: Select best matches and write'm out
        # TODO: Align the matches
        # TODO: Un-translate the sequences
        # TODO: Pack the results and clean the house
