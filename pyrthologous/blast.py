"""Wrapper for blast program."""

import os
import subprocess
from .config import BLASTP, BLAST_DB_MAKER, OUTPUT


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
               "-out", tgt,
               "-dbtype", "prot",
               "-logfile", "DBBLAST.log"]

    p = subprocess.Popen(command)
    p.communicate()

    return tgt


def reciprocal_blastp(query_subject):
    """Blast the first of query_subject against the second and viceversa."""

    for query, subject in query_subject, query_subject[::-1]:
        # Make the database
        db = make_blast_db(subject,
                           os.path.join(os.path.dirname(subject), OUTPUT))
        # Blast'em
        stdout, stderr = blastp(query, subject)
        # TODO: delete the database
        # Yield and repeat
        yield stdout, stderr
