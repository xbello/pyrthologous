"""Run all the functions with the genomes."""

import os
import subprocess
from Bio import SeqIO
from tempfile import NamedTemporaryFile

import blast
from .config import *
from .prepare import translate_fasta
from .utils import detranslate


def align(pair):
    """Return a string with the file name to the alignment between a pair."""
    # Save the matches to temporary files
    tmp_align = NamedTemporaryFile()
    align_file = tmp_align.name + ".al.fasta"

    SeqIO.write(pair.values(), tmp_align, "fasta")
    # Rewind the file so Muscle can read it
    tmp_align.seek(0)
    # Align the matches
    proc = subprocess.Popen(
        [MUSCLE,
         "-in", tmp_align.name,
         "-out", align_file],
        stderr=subprocess.PIPE)

    _, err = proc.communicate()
    tmp_align.close()

    # Check if Muscle failed
    if "ERROR" in err:
        raise IOError(err)

    return align_file


def blast_pair(pair):
    """Return a dict with the results of a reciprocal blast."""

    out_dicts = []
    abs_paths = [os.path.join(BASE_PATH, OUTPUT, x) for x in pair]

    outputs, errs = zip(*blast.reciprocal_blastp(abs_paths))
    # Put each output in a list
    for output in outputs:
        out_dict = {}
        for k, v in [blast.blastout_to_tuple(x) for x in output]:
            out_dict[k] = v
        out_dicts.append(out_dict)

    return out_dicts


def get_seq_from(pair, dict_of_seqs):
    """Return a dict with the sequence for the key name."""

    pair_dict = dict.fromkeys(pair)
    for seq in pair:
        pair_dict[seq] = dict_of_seqs[0].get(seq) or dict_of_seqs[1].get(seq)

    return pair_dict


def get_best_matches(dict_of_matches):
    """Return a list with a pair of reciprocal matches."""
    matches = []

    for k, v in dict_of_matches[0].iteritems():
        if v[0] in dict_of_matches[1] and k == dict_of_matches[1][v[0]][0]:
            matches.append((k, v[0]))

    return matches


def init_pair(pair):
    """Return the genomes SeqIO'ed as a dict from the bases and aa."""

    genomes = {"base": [],
               "aa": []}

    for specie in pair:
        aa_genome = os.path.join(BASE_PATH, OUTPUT, specie)
        base_genome = os.path.join(BASE_PATH, GENOMES, specie)

        if not os.path.isfile(base_genome):
            raise IOError("File {0} doesn't exist".format(specie))

        # Translate both genomes only if they doesn't exist
        if not os.path.isfile(aa_genome):
            translate_fasta(base_genome,
                            output_path=os.path.dirname(aa_genome))
        # If we are here without error, append each genome to genomes
        genomes["base"].append(
            SeqIO.to_dict(SeqIO.parse(base_genome, "fasta")))
        genomes["aa"].append(
            SeqIO.to_dict(SeqIO.parse(aa_genome, "fasta")))

    return genomes


if __name__ == "__main__":
    for pair in COMPARE:
        # Generate the genomes
        import sys
        genomes = init_pair(pair)

        best_matches = get_best_matches(blast_pair(pair))

        # Align each match
        alignments = []
        for match in best_matches:
            align_file = align(get_seq_from(match, genomes["aa"]))
            # De-translate each alignment
            base_align = detranslate(
                align_file, get_seq_from(match, genomes["base"]))
            SeqIO.write(base_align, open(align_file, "w"), "fasta")
            print align_file
            sys.exit()


        # TODO: Pack the results and clean the house
