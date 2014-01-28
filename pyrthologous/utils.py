"""Miscelaneous functions for the module."""

import itertools
import os
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord


def translate(seqs):
    """Return the translated sequence in seqs Seq records."""

    for sequence in seqs:
        new_record = SeqRecord(sequence.seq.translate(),
                               id=sequence.id,
                               description="")
        yield new_record
