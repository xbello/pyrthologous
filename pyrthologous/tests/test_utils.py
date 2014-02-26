"""Test the misc functions in utils."""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from unittest import TestCase
from .. import utils


class testUtils(TestCase):

    def setUp(self):
        self.seq1 = SeqRecord(Seq("ATG"), id="Single Amino")
        self.seq2 = SeqRecord(Seq("AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"),
                              id="Multiple Aminos")
        self.seq3 = SeqRecord(Seq("AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"),
                              id="Bad header,")

    def test_translate(self):

        recs = [x for x in utils.translate([self.seq1, self.seq2])]

        self.assertEqual(str(recs[0].seq), "M")
        self.assertEqual(str(recs[1].seq), "KEREREQERCSS")

        self.assertEqual(recs[0].id, "Single Amino")

    def test_reverse_translate(self):
        seq_bases = "AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"
        seq_aa = "K-ERERE-QERCSS"

        self.assertEqual(
            utils.reverse_translate(seq_aa, seq_bases),
            "AAG---GAGAGGGAGAGAGAG---CAGGAACGATGTTCTTCC")

    def test_clean_title(self):
        good = "This_is_a_good_title"
        good1 = good + " with some additions, and other"
        bad = good + ", or nasty commas"

        self.assertEqual(utils.clean_title(good), (good, "", good))
        self.assertEqual(utils.clean_title(good1), (good, "", good1))
        self.assertEqual(utils.clean_title(bad), (good, "", bad))
