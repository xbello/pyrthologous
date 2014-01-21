"""Test for prepare.py module."""

import os
from collections import Counter
from unittest import TestCase

from .. import prepare


class testPrepare(TestCase):
    def setUp(self):
        self.genome = "mock_genome"
        self.genome_path = os.path.join(os.path.dirname(__file__), self.genome)
        self.fastas = ["1.fas", "2.fas", "3.fas",
                       "prot1.fasta", "prot2.fasta"]

    def test_check_compare(self):
        self.assertTrue(prepare.check_compare(path=self.genome_path,
                                              compare=(("1.fas", "2.fas"),
                                                       ("2.fas", "3.fas"))
                                              ))
        self.assertRaises(AssertionError,
                          prepare.check_compare,
                          *[],
                          **{"path": self.genome_path,
                             "compare": (("4.fas", "2.fas"),)})

    def test_get_fasta_names(self):

        self.assertEqual(Counter(prepare.get_fasta_names(self.genome_path)),
                         Counter(self.fastas))

    def test_translate_fasta(self):
        genome_path = os.path.join(os.path.dirname(__file__),
                                   "mock_genome")

        fasta_path = os.path.join(genome_path, "1.fas")
        output_file = os.path.join(genome_path, "output", "1.fas")

        if not os.path.exists(os.path.dirname(output_file)):
            os.mkdir(os.path.dirname(output_file))

        self.assertFalse(os.path.isfile(output_file))
        prepare.translate_fasta(fasta_path,
                                os.path.join(genome_path, "output"))
        self.assertTrue(os.path.isfile(output_file))

        trans_file = open(output_file, "r").readlines()

        self.assertEqual(trans_file[0].strip(),
                         ">6_BEWA_012860_804_signal")
        self.assertEqual(trans_file[1].strip()[:20], "MKKCVSSSIIALFTLFISTG")
        self.assertEqual(trans_file[-1].strip()[:20], "AYNEARLVSKHINKLLEVDV")

        os.unlink(output_file)
