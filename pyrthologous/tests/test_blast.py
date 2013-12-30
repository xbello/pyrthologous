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

        self.blast_output = [
            "".join(["Protein_query\tProt2\t99.33\t300\t0\t2",
                     "\t122\t420\t2\t300\t0.0\t 626"]),
            "".join(["Protein_query\tProt2\t21.88\t64\t30\t1",
                     "\t74\t137\t175\t218\t1.8\t16.9"]),
            "".join(["Protein_query\tProt2\t33.33\t18\t12\t0",
                     "\t37\t54\t331\t348\t6.9\t15.4"])]

    def tearDown(self):
        for filep in os.listdir(self.tgt_path):
            if filep.endswith((".phr", ".pin", ".psq")):
                os.unlink(os.path.join(self.tgt_path, filep))

    def test_blastp(self):
        db = blast.make_blast_db(self.prots,
			os.path.join(self.tgt_path, "prots"))

	blast_stdout, blast_stderr = blast.blastp(self.seq1, db)

        self.assertEqual(blast_stdout.splitlines(), self.blast_output)
	self.assertEqual(blast_stderr, "")

    def test_make_blast_db(self):
        self.assertEqual(
            blast.make_blast_db(self.prots, os.path.join(self.tgt_path,
                                                         "prots")),
            os.path.join(self.tgt_path, "prots"))
