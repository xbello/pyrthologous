"""Miscelaneous functions for the module."""

import logging
import subprocess
from tempfile import NamedTemporaryFile

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def align(pair, config):
    """Return a string with the file name to the alignment between a pair."""
    # Save the matches to temporary files
    tmp_align = NamedTemporaryFile()
    align_file = tmp_align.name + ".al.fasta"

    SeqIO.write(pair.values(), tmp_align, "fasta")
    # Rewind the file so Muscle can read it
    tmp_align.seek(0)
    # Align the matches
    logger.info("Aligning {0} vs {1}.".format(*pair.keys()))

    proc = subprocess.Popen(
        [config.MUSCLE,
         "-in", tmp_align.name,
         "-out", align_file],
        stderr=subprocess.PIPE)

    _, err = proc.communicate()
    tmp_align.close()

    # Check if Muscle failed
    if "ERROR" in err:
        logger.critical("Failed while aligning {0}".format(err))
        raise IOError(err)

    return align_file


def clean_dict(SeqIO_dict):
    """Return a dict with the keys cleaned."""

    for d in SeqIO_dict:
        for k in d.keys():
            if k.endswith(","):
                # Delete awful last commas
                d[k[:-1]] = d[k]
                d[k[:-1]].id = d[k[:-1]].id[:-1]
                del d[k]

    return SeqIO_dict


def detranslate(filename, base_seq):
    """Return the de-translated sequences of filename."""

    alignment = SeqIO.to_dict(AlignIO.read(filename, "fasta"))
    base_alignment = []

    for k, v in alignment.items():
        base_alignment.append(
            SeqRecord(reverse_translate(v.seq, base_seq[k].seq),
                      id=k,
                      description=""))

    return base_alignment


def get_seq_from(pair, dict_of_seqs):
    """Return a set with the sequences for the keys names in pair."""

    # XXX If both sequences in pair are the same, dict will have only one key
    # XXX Should change this from dicts to sets of two sequences.
    pair_dict = dict.fromkeys(pair)
    for seq in pair:
        pair_dict[seq] = dict_of_seqs[0].get(seq) or dict_of_seqs[1].get(seq)
        if not pair_dict[seq]:
            raise KeyError("{0} not found in genomes.".format(seq))
        if not pair_dict[seq].id:
            raise KeyError("{0} didn't have an id field.".format(seq))

    return pair_dict


def reverse_translate(aa_seq, base_seq):
    """Return a string mapped from seq_aa to bases."""

    new_seq = ""
    from_p = 0

    for pos in aa_seq:
        if pos == "-":
            new_seq += "---"
        else:
            new_seq += base_seq[from_p:from_p + 3]
            from_p += 3

    return new_seq


def translate(seqs):
    """Return the translated sequence in seqs Seq records."""

    for sequence in seqs:
        new_record = SeqRecord(sequence.seq.translate(),
                               id=sequence.id,
                               description="")
        yield new_record
