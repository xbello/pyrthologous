"""Test the misc functions in utils."""

from unittest import TestCase
from collections import Counter
from os.path import abspath, join
from .. import utils

import test_config as tc


class testUtils(TestCase):

    def test_make_reciprocals(self):
        self.assertEqual(
            utils.make_reciprocals(["1.fas", "2.fas", "3.fas"]),
            [('1.fas', '2.fas'), ('1.fas', '3.fas'), ('2.fas', '3.fas')])

    def test_normalize(self):
        normalized = [x for x in utils.normalize("tests/mock_genome/1.fas")]

        self.assertEqual(
            normalized[0],
            (1, "6_BEWA_012860_804_signal"))

        self.assertEqual(5, len(normalized))

    def test_reciprocal_genomes(self):
        self.assertEqual(
            Counter(
                [('1.fas', '2.fas'), ('1.fas', '3.fas'), ('2.fas', '3.fas')]),
            Counter([tuple(sorted(x)) for x in
                     utils.reciprocal_genomes("tests/mock_genome", "fas")]))

    def test_select_genomes(self):
        self.assertEqual(
            Counter(utils.select_genomes(tc.BASE_PATH, "fas")),
            Counter(["1.fas", "2.fas", "3.fas"]))
        self.assertFalse(
            utils.select_genomes("tests/this_path_is_fake", "fas"))

    def test_tag_fasta(self):
        self.assertEqual(
            utils.tag_fasta(join(tc.BASE_PATH, "1.fas")),
            abspath("tests/mock_genome/1_{0}.fas".format(tc.TAG)))
        self.assertEqual(
            utils.tag_fasta(join(tc.BASE_PATH, "1.fas"), tag="MOCK"),
            abspath("tests/mock_genome/1_MOCK.fas"))

    def test_translate(self):
        self.assertEqual(
            utils.translate("ATG"), "M")
        self.assertEqual(
            utils.translate("AAGGAGAGGGAGAGAGAGCAGGAACGATGTTCTTCC"),
            "KEREREQERCSS")
