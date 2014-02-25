"""Run all the functions with the genomes."""

import os
import tarfile
from Bio import SeqIO

from pyrthologous import blast
from pyrthologous.prepare import translate_fasta
from pyrthologous.utils import align, detranslate, clean_dict


def blast_pair(pair):
    """Return a dict with the results of a reciprocal blast."""

    out_dicts = []
    abs_paths = [os.path.join(c.BASE_PATH, c.OUTPUT, x) for x in pair]

    outputs, errs = zip(*blast.reciprocal_blastp(abs_paths, c))

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
        if not pair_dict[seq]:
            raise KeyError("{0} not found in genomes.".format(seq))
        if not pair_dict[seq].id:
            raise KeyError("{0} didn't have an id field.".format(seq))

    return pair_dict


def get_best_matches(dict_of_matches):
    """Yield a list with a pair of reciprocal matches."""

    for k, v in dict_of_matches[0].iteritems():
        if v[0] in dict_of_matches[1] and k == dict_of_matches[1][v[0]][0]:
            yield (k, v[0])


def init_pair(pair):
    """Return the genomes SeqIO'ed as a dict from the bases and aa."""

    genomes = {"base": [],
               "aa": []}

    for specie in pair:
        aa_genome = os.path.join(c.BASE_PATH, c.OUTPUT, specie)
        base_genome = os.path.join(c.BASE_PATH, c.GENOMES, specie)

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

        genomes["base"] = clean_dict(genomes["base"])
        genomes["aa"] = clean_dict(genomes["aa"])

    return genomes


def main(pair, output_path):
    """Run the main script, useful for multicoring."""

    # Generate the pair id
    pair_name = "{0}_vs_{1}".format(
        *[x.split(".")[0] for x in pair])

    # Generate the genomes
    genomes = init_pair(pair)

    best_matches = get_best_matches(blast_pair(pair))

    # Create the tar pack
    tar_file = tarfile.open(
        os.path.join(output_path, pair_name + ".tgz"), "w:gz")

    # Align each match
    for match in best_matches:
        align_file = align(get_seq_from(match, genomes["aa"]), c)
        # De-translate each alignment
        base_align = detranslate(
            align_file, get_seq_from(match, genomes["base"]))
        SeqIO.write(base_align, open(align_file, "w"), "fasta")
        # Add each alignment to the tar pack
        tar_file.add(align_file, recursive=False)
        # Remove the temporary file
        os.unlink(align_file)

    tar_file.close()
    print tar_file.name

if __name__ == "__main__":
    # Get the config for a config.py of CWD dir
    if os.path.isfile(os.path.join(os.getcwd(), "config.py")):
        # There is a config.py file in the calling directory, load it
        from imp import load_source
        c = load_source("config",
                        os.path.join(os.getcwd(), "config.py"))
        # Python3:
        # from importlib.machinery import SourceFileLoader
        # c = SourceFileLoader("config", "config.py").load_module()

    else:
        import pyrthologous.config as c

    from multiprocessing import Pool
    pool = Pool(4)
    threads = []

    for pair in c.COMPARE:
        output_path = os.path.join(c.BASE_PATH, c.OUTPUT)
        threads.append(
            pool.apply_async(main, args=[pair, output_path]))

    pool.close()
    pool.join()
