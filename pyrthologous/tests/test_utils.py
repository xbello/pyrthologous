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

    def test_clean_dict(self):
        list_seqs = [{"Bad header,": self.seq3},
                     {"Multiple Aminos": self.seq2}]

        clean_dict = utils.clean_dict(list_seqs)
        # Keys has been cleaned
        self.assertTrue("Bad header" in clean_dict[0].keys())
        # Ids has been cleaned
        self.assertEqual(
            clean_dict[0]["Bad header"].id,
            "Bad header")
