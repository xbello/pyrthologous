"""Test the misc functions in utils."""

from unittest import TestCase
from collections import Counter
from .. import utils


class testUtils(TestCase):

    def test_make_reciprocals(self):
        self.assertEqual(
            utils.make_reciprocals(["1.fas", "2.fas", "3.fas"]),
            [('1.fas', '2.fas'), ('1.fas', '3.fas'), ('2.fas', '3.fas')])

    def test_reciprocal_genomes(self):
        self.assertEqual(
            Counter(
                [('1.fas', '2.fas'), ('1.fas', '3.fas'), ('2.fas', '3.fas')]),
            Counter([tuple(sorted(x)) for x in
                     utils.reciprocal_genomes("tests/mock_genome", "fas")]))

    def test_select_genomes(self):
        self.assertEqual(
            Counter(utils.select_genomes("tests/mock_genome", "fas")),
            Counter(["1.fas", "2.fas", "3.fas"]))
        self.assertFalse(
            utils.select_genomes("tests/this_path_is_fake", "fas"))
