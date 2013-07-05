#!/usr/bin/env python
import os
from Bio import SeqIO

from config import SEP, MUSCLE

muscle = "/home/xbello/bin/muscle"

def comp_dicts(dict_1, dict_2):
    '''Dict, Dict -> Dict

    If two dicts contains reciprocal keys, add that values to a new dict'''

    d = {}

    for k, v in dict_1.iteritems():
        if dict_2[v] == k:
            d[k] = v

    return d

def get_dir(tab_blast_output):
    '''String -> Dict {String: String}

    From a tabular blast makes a dict with query: target'''

    d = {}

    with open(tab_blast_output, "r") as tab_blast:
        for line in tab_blast:
            data = line.split("\t")
            d[data[0]] = data[1]
            
    return d

def get_seq_from_fasta(seq_id, fasta_file):
    '''String -> String

    Given a seq id, extract the sequence in FASTA file.'''

    records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    return records[seq_id]

def make_align(seq_id_one, seq_id_two):
    '''String, String -> String

    Given two ids, compute an alignment (MUSCLE) and return the path of that
    alignment'''

    seq_one = get_seq_from_fasta(seq_id_one, fasta_src_1)
    seq_two = get_seq_from_fasta(seq_id_two, fasta_src_2)
    output_file = os.path.join(os.getcwd(), align_dir, "{0}_{1}.fas".format(
        seq_one.name, seq_two.name))
    handle = open(output_file, "w")
    SeqIO.write([seq_one, seq_two], handle, "fasta")
    handle.close()

    os.system("{0} -in {1} -out {2}".format(
        MUSCLE, output_file, output_file + ".al"))

    return output_file + ".al"

def make_transtable(transtable_filename):
    '''String -> Dict

    Given a transtable with the format:

      Long name of sequence,number

    Return a dict with the format:

      {"number": "Long name of sequence"}
    '''
    d = {}

    with open(transtable_filename, "r") as t_f:
        for line in t_f:
            data = line.strip().split(SEP)
            #Some lines got commas in the description
            d[data[-1]] = "".join(data[0:-1])

    return d

def protein_to_bases(alignment, seq_id1, seq_id2):
    '''String -> Bool

    Given the name of a file, de-translate the sequences from aa to bases'''

    records_AA = [x for x in SeqIO.parse(alignment, "fasta")]

    records_BASES_1 = SeqIO.to_dict(SeqIO.parse(bases_src_1, "fasta"))
    records_BASES_2 = SeqIO.to_dict(SeqIO.parse(bases_src_2, "fasta"))
    output_file = os.path.join(os.getcwd(), align_dir,
        "{0}_{1}_bases.fas".format(seq_id1, seq_id2))

    handle = open(output_file, "w")

    seq_1 = records_AA[0]
    seq_2 = records_AA[1]
    new_record_name_1 = trans_table_1[seq_id1]
    new_record_name_2 = trans_table_2[seq_id2]

    for name, seq in records_BASES_1.iteritems():
        if name.startswith(new_record_name_1):
            handle.write(">{0}\n".format(name))
            handle.write("{0}\n".format(seq_aa_to_bases(seq_1.seq, seq.seq)))

    for name, seq in records_BASES_2.iteritems():
        if name.startswith(new_record_name_2):
            handle.write(">{0}\n".format(name))
            handle.write("{0}\n".format(seq_aa_to_bases(seq_2.seq, seq.seq)))

    handle.close()

    return True
    
def seq_aa_to_bases(seq_aa, seq_bases):
    '''String, String -> String

    Given a sequence in aa and "-", substitute each aa with a triplet'''

    new_seq = ""
    from_p = 0
    to_p = 3

    for pos in seq_aa:
        if pos == "-":
            new_seq += "---"
        else:
            new_seq += seq_bases[from_p:to_p]
            from_p += 3
            to_p += 3

    return new_seq

if __name__ == "__main__":
    import config
    import utils
    import sys

    for x in config.LEAF_PATHS:
        for pair in utils.make_reciprocals(x):
            organism_1 = "_".join(pair[0].split(".")[0].split("_")[:2])
            organism_2 = "_".join(pair[1].split(".")[0].split("_")[:2])

            genus = organism_1.split("_")[0]

            organism_1_abbr = "".join([x[0] for x in organism_1.split("_")])
            organism_2_abbr = "".join([x[0] for x in organism_2.split("_")])
        
            file_1 = os.path.join(
                config.OUTPUT,
                "_vs_".join([organism_1_abbr, organism_2_abbr]) + ".blast")
            file_2 = os.path.join(
                config.OUTPUT,
                "_vs_".join([organism_2_abbr, organism_1_abbr]) + ".blast")

            fasta_src_1 = os.path.join(
                config.BASE_PATH,
                genus,
                pair[0])
            bases_src_1 = fasta_src_1.replace("_NORM", "") + "ta"
            fasta_src_2 = os.path.join(
                config.BASE_PATH,
                genus,
                pair[1])
            bases_src_2 = fasta_src_2.replace("_NORM", "") + "ta"

            align_dir = "_vs_".join([organism_1_abbr, organism_2_abbr])
            if not os.path.isdir(align_dir): os.mkdir(align_dir)
           
            trans_table_1 = make_transtable(bases_src_1.replace("fasta", "txt"))
            trans_table_2 = make_transtable(bases_src_2.replace("fasta", "txt"))

            one_to_two = get_dir(file_1)
            two_to_one = get_dir(file_2)

            new_dict = comp_dicts(one_to_two, two_to_one)

            for k, v in new_dict.iteritems():
                alignment_name = make_align(k, v)
                protein_to_bases(alignment_name, k, v)
