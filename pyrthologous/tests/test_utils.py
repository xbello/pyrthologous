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
