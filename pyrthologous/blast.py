"""Wrapper for blast program."""

import os
import subprocess
from pyrthologous.config import BLASTP, BLAST_DB_MAKER


def blastp(query, subject):
    """Return the output for a blast from query to subject."""

    command = [BLASTP,
               "-query", query,
               "-db", subject,
               "-outfmt", "6",
               "-max_target_seqs", "1"]

    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    return stdout, stderr


def make_blast_db(src, tgt):
    """Return the Database name created from ``src``."""

    command = [BLAST_DB_MAKER,
               "-in", src,
               "-out", os.path.join(os.getcwd(), tgt),
               "-dbtype", "prot",
               "-logfile", "DBBLAST.log"]

    p = subprocess.Popen(command)
    p.communicate()

    return tgt
