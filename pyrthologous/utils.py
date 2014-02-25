"""Miscelaneous functions for the module."""

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord


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
