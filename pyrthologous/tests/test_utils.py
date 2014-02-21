"""Test the misc functions in utils."""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from unittest import TestCase
from .. import utils


class testUtils(TestCase):

    def test_translate(self):
        seq1 = SeqRecord(Seq("ATG"), id="Single Amino")
        seq2 = SeqRecord(Seq("AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"),
                         id="Multiple Aminos")

        recs = [x for x in utils.translate([seq1, seq2])]

        self.assertEqual(str(recs[0].seq), "M")
        self.assertEqual(str(recs[1].seq), "KEREREQERCSS")

        self.assertEqual(recs[0].id, "Single Amino")

    def test_reverse_translate(self):
        seq_bases = "AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"
        seq_aa = "K-ERERE-QERCSS"

        self.assertEqual(
            utils.reverse_translate(seq_aa, seq_bases),
            "AAG---GAGAGGGAGAGAGAG---CAGGAACGATGTTCTTCC")
