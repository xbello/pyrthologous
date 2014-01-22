"""Test for module blast."""

import os
from unittest import TestCase

from .. import blast


class testBlast(TestCase):
    def setUp(self):
        self.tgt_path = "tests/mock_genome"
        self.prots = os.path.abspath(
            os.path.join(self.tgt_path, "prot1.fasta"))
        self.seq1 = os.path.abspath(
            os.path.join(self.tgt_path, "prot2.fasta"))

        self.blast_output = "Protein_query\tProt2\t99.33\t300\t0\t2" +\
            "\t122\t420\t2\t300\t0.0\t 626"

    def tearDown(self):
        for filep in os.listdir(self.tgt_path):
            if filep.endswith((".phr", ".pin", ".psq")):
                os.unlink(os.path.join(self.tgt_path, filep))

    def test_blastp(self):
        """Test blastp output: shown blastp output vs expected output."""
        db = blast.make_blast_db(self.prots,
                                 os.path.join(self.tgt_path, "prots"))

        stdout, stderr = blast.blastp(self.seq1, db)

        self.assertEqual(stdout.split("\n")[0], self.blast_output)
        self.assertEqual(stderr, "")

    def test_make_blast_db(self):
        self.assertEqual(
            blast.make_blast_db(self.prots, os.path.join(self.tgt_path,
                                                         "prots")),
            os.path.join(self.tgt_path, "prots"))
