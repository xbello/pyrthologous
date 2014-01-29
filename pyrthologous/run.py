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
            this_genome = os.path.join(BASE_PATH, OUTPUT, specie)
            if not os.path.isfile(this_genome):
                translate_fasta(os.path.join(BASE_PATH, GENOMES, specie),
                                output_path=os.path.dirname(this_genome))

        # Reciprocal blast'em
        abs_paths = [os.path.join(BASE_PATH, OUTPUT, x) for x in pair]

        outputs, errs = zip(*blast.reciprocal_blastp(abs_paths))
        # Put each output in a list
        out_dicts = []
        for output in outputs:
            out_dict = {}
            for k, v in [blast.blastout_to_tuple(x) for x in output]:
                out_dict[k] = v
            out_dicts.append(out_dict)

        # TODO: Select best matches and write'm out
        matches = []
        for k, v in out_dicts[0].items():
            if v[0] in out_dicts[1] and k == out_dicts[1][v[0]][0]:
                matches.append([k, v[0]])
                # TODO: Do the backpass

        # TODO: This prints BUT should extract and write to file
        print matches



        # TODO: Align the matches
        # TODO: Un-translate the sequences
        # TODO: Pack the results and clean the house
