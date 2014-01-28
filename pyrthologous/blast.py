"""Wrapper for blast program."""

import __builtin__
import os
import subprocess
from operator import itemgetter
from .config import BLASTP, BLAST_DB_MAKER, OUTPUT

# Globals for column positions of tab-blast output
QUERY = 0
IDENTITY = 2


def blastp(query, subject):
    """Return the output for a blast from query to subject."""

    command = [BLASTP,
               "-query", query,
               "-db", subject,
               "-outfmt", "6"]
               # XXX This flag have erratical behaviour among BLAST versions
               # 2.2.27+ output differs from 2.2.28+
               #"-max_target_seqs", "1"]

    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if stderr:
        raise IOError(stderr)

    stdout = get_best_from_blast_output(stdout)

    return stdout, stderr


def make_blast_db(src, tgt):
    """Return the Database name created from ``src``."""

    if not os.path.isdir(tgt):
        os.mkdir(tgt)
    tgt = os.path.join(tgt, os.path.basename(src))

    command = [BLAST_DB_MAKER,
               "-in", src,
               "-out", tgt,
               "-dbtype", "prot",
               "-logfile", "DBBLAST.log"]

    p = subprocess.Popen(command)
    p.communicate()

    return tgt


def reciprocal_blastp(query_subject):
    """Blast the first of query_subject against the second and viceversa.

    query_subject is a set as (query, subject)

    """

    for query, subject in query_subject, query_subject[::-1]:
        # Make the database
        db = make_blast_db(subject,
                           os.path.join(os.path.dirname(subject), OUTPUT))

        # Blast'em
        stdout, stderr = blastp(query, db)

        # Yield and repeat
        yield stdout, stderr


def listify_blast_output(blast_output, casts=[]):
    """Return a list wit the blast output."""

    return_list = []

    for match_line in blast_output.rstrip().split("\n"):
        new_line = match_line.split("\t")
        for cast in casts:
            # Cast the columns in "casts" to types for further comparisons
            new_line[cast[0]] = getattr(
                __builtin__, cast[1])(new_line[cast[0]])
        return_list.append(new_line)

    return return_list


def get_best_from_blast_output(blast_output):
    """Return only the best matches for each group of matches."""

    blast_list = listify_blast_output(
        blast_output, casts=[(IDENTITY, "float")])

    return "\n".join(
        ["\t".join([str(x) for x in y])
         for y in simplify_blast_output(blast_list=blast_list, group=[])])


def get_best_from_group(blast_list, position):
    """Return the best item from a list."""

    # Sort the matches by position and get the best of them
    best_match = sorted(
        blast_list, key=itemgetter(position), reverse=True)[0]

    return best_match


def simplify_blast_output(blast_list=[], group=[]):
    """Generator that yield the best result for each group of lines."""

    # If there are lines left to process...
    while blast_list:
        # Pop the first line in the batch
        this_line = blast_list.pop(0)
        group.append(this_line)

        if not blast_list or (blast_list[0][QUERY] != this_line[QUERY]):
            # List has been just exhausted or the next line is a new group
            # Either two options: process current group, yield and continue
            yield get_best_from_group(group, IDENTITY)
            group = []

        # Recurse the function while blast_list isn't exhausted
        for x in simplify_blast_output(blast_list=blast_list, group=group):
            yield x
