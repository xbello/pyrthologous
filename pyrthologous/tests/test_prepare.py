import os
from collections import Counter
from unittest import TestCase

from .. import prepare

class testPrepare(TestCase):
    def setUp(self):
        self.genome = "mock_genome"
        self.fastas = ["1.fas", "2.fas", "3.fas",
                       "prot1.fasta", "prot2.fasta"]

    def test_get_fasta_names(self):
        path = os.path.join(os.path.dirname(__file__), self.genome)

        self.assertEqual(Counter(prepare.get_fasta_names(path)),
                         Counter(self.fastas))
