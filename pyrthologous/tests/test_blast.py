"""Test for module blast."""

import os
import shutil
from imp import load_source
from unittest import TestCase

from .. import blast
from .. import config


class testBlast(TestCase):
    def setUp(self):
        self.tgt_path = os.path.join(os.path.dirname(__file__), "mock_genome")
        # This is a fasta with two protein sequences, named Prot1 and Prot2
        self.prots = os.path.abspath(
            os.path.join(self.tgt_path, "prot1.fasta"))
        # This is a fasta with only one protein sequence, Protein_query
        self.prot2 = os.path.abspath(
            os.path.join(self.tgt_path, "prot2.fasta"))

        self.config = load_source("config", "config.py")

        self.blast_output = "Protein_query\tProt2\t99.33\t300\t0\t2" +\
            "\t122\t420\t2\t300\t0.0\t 626"
        self.blast_output_two_best = \
            "Prot1\tProtein_query\t31.25\t16\t11\t0\t168\t183\t262\t" +\
            "277\t3.3\t15.0\n" +\
            "Prot2\tProtein_query\t99.33\t300\t0\t2\t2\t300\t122\t" +\
            "420\t0.0\t 626"

        self.blast_output2 = "Protein_query\tProt2\t99.33\t300\t0\t2" +\
            "\t122\t420\t2\t300\t0.0\t 626\n" +\
            "Protein_query\tProt2\t21.88\t64\t30\t1\t74\t137\t175\t218" +\
            "\t1.8\t16.9\n" +\
            "Protein_query\tProt2\t33.33\t18\t12\t0\t37\t54\t331\t348" +\
            "\t6.9\t15.4"

        self.blast_output3 = "Prot1\tProtein_query\t27.78\t54\t29\t2\t" +\
            "137\t185\t359\t4070.35\t18.1\n" +\
            "Prot1\tProtein_query\t31.25\t16\t11\t0\t168\t183\t262\t" +\
            "277\t3.3\t15.0\n" +\
            "Prot1\tProtein_query\t25.58\t43\t32\t0\t34\t76\t238\t" +\
            "280\t4.9\t14.6\n" +\
            "Prot2\tProtein_query\t99.33\t300\t0\t2\t2\t300\t122\t" +\
            "420\t0.0\t 626\n" +\
            "Prot2\tProtein_query\t21.88\t64\t30\t1\t175\t218\t74\t" +\
            "137\t1.1\t16.9\n" +\
            "Prot2\tProtein_query\t33.33\t18\t12\t0\t331\t348\t37\t" +\
            "54\t4.0\t15.4"

    def tearDown(self):
        for filep in os.listdir(self.tgt_path):
            if filep.endswith((".phr", ".pin", ".psq")):
                os.unlink(os.path.join(self.tgt_path, filep))
            if os.path.isdir(os.path.join(self.tgt_path, filep)):
                if filep in ["prots", config.OUTPUT]:
                    shutil.rmtree(os.path.join(self.tgt_path, filep))

    def test_blastout_to_tuple(self):
        blast_dict = {}
        for k, v in [blast.blastout_to_tuple(x.split("\t")) for x in
                     self.blast_output_two_best.split("\n")]:
            blast_dict[k] = v

        self.assertEqual(
            blast_dict,
            {"Prot1": ("Protein_query",
                       ["31.25", "16", "11", "0", "168", "183",
                        "262", "277", "3.3", "15.0"]),
             "Prot2": ("Protein_query",
                       ["99.33", "300", "0", "2", "2", "300",
                        "122", "420", "0.0", " 626"])})

    def test_blastp(self):
        """Test blastp output: shown blastp output vs expected output."""
        db = blast.make_blast_db(self.prots,
                                 os.path.join(self.tgt_path, "prots"),
                                 self.config)

        stdout, stderr = blast.blastp(self.prot2, db, self.config)

        self.assertEqual("\t".join([x for x in stdout][0]),
                         self.blast_output)
        self.assertEqual(stderr, "")

    def test_make_blast_db(self):
        self.assertEqual(
            blast.make_blast_db(self.prots,
                                os.path.join(self.tgt_path, "prots"),
                                self.config),
            os.path.join(self.tgt_path,
                         "prots",
                         os.path.basename(self.prots)))

    def test_reciprocal_blastp_outputs(self):
        pair = (self.prots, self.prot2)
        stdouts, stderrs = zip(*blast.reciprocal_blastp(pair, self.config))

        # Assert no errors yield by blastp
        self.assertEqual(stderrs, ("", ""))

        # Assert the output are correct
        stdout_1 = ["\t".join(x) for x in stdouts[0]]
        stdout_2 = ["\t".join(x) for x in stdouts[1]]
        self.assertItemsEqual(
            [stdout_1, stdout_2],
            [self.blast_output.split("\n"),
             self.blast_output_two_best.split("\n")])

    def test_listify_blast_line(self):
        # Only one line, no casting
        list_blast_output = self.blast_output.split("\t")
        self.assertEqual(
            blast.listify_blast_line(self.blast_output),
            list_blast_output)

        # Cast the second column to a float
        list_blast_output[2] = float(list_blast_output[2])
        self.assertEqual(
            blast.listify_blast_line(self.blast_output, casts=[(2, "float")]),
            list_blast_output)

    def test_get_best_from_blast_output(self):
        # Just one group
        self.assertEqual(
            [x for x in
             blast.get_best_from_blast_output(self.blast_output2)][0],
            self.blast_output.split("\t"))

        #Multiple groups
        self.assertItemsEqual(
            [x for x in blast.get_best_from_blast_output(
                self.blast_output3)],
            [x.split("\t") for x in self.blast_output_two_best.split("\n")])

    def test_simplify_blast_output(self):
        # One group get simplified to its best line
        output_as_list = [x.split("\t") for x in
                          self.blast_output2.split("\n")]

        bests = [x for x in blast.simplify_blast_output(
            blast_list=output_as_list, group=[])]

        self.assertEqual(bests, [self.blast_output.split("\t")])

        # Multiple groups simplified to its best lines
        output_as_list = [x.split("\t") for x in
                          self.blast_output3.split("\n")]

        bests = [x for x in blast.simplify_blast_output(
            blast_list=output_as_list, group=[])]

        self.assertEqual(
            bests,
            [x.split("\t") for x in self.blast_output_two_best.split("\n")])
