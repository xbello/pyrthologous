"""Wrapper for blast program."""

import __builtin__
import logging
import os
import subprocess
from operator import itemgetter

# Globals for column positions of tab-blast output
QUERY = 0
IDENTITY = 2

logger = logging.getLogger(__name__)


def blastout_to_tuple(blast_line):
    """Return a tuple (key, value) with each blast_line entry."""

    #data = blast_line.rstrip().split("\t")
    return blast_line[0], (blast_line[1], blast_line[2:])


def blastp(query, subject, config):
    """Return the output for a blast from query to subject."""

    command = [config.BLASTP,
               "-query", query,
               "-db", subject,
               "-outfmt", "6",
               # XXX -max_target_seqs flag have erratical behaviour among
               # BLAST versions. 2.2.27+ output differs from 2.2.28+
               "-max_target_seqs", "1"]

    logger.info("Trying to run blastp for {0} against {1}.".format(
        subject, query))
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if stderr:
        logger.critical("Failed blastp with '{0}'.".format(stderr))
        raise IOError(stderr)

    logger.info("Ended blastp.")

    stdout = get_best_from_blast_output(stdout)

    return stdout, stderr


def make_blast_db(src, tgt, config):
    """Return the Database name created from ``src``."""

    if not os.path.isdir(tgt):
        os.mkdir(tgt)
    tgt = os.path.join(tgt, os.path.basename(src))

    command = [config.BLAST_DB_MAKER,
               "-in", src,
               "-out", tgt,
               "-dbtype", "prot"]

    logger.info("Compiling database with {0}.".format(src))
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if stderr:
        logger.critical(stderr)

    return tgt


def reciprocal_blastp(query_subject, config):
    """Return a list of tuples (output generator, errors).

    Blast the first of query_subject against the second and viceversa.

    query_subject is a set as (query, subject)

    """

    outputs = []
    for query, subject in query_subject, query_subject[::-1]:
        # Make the database
        db = make_blast_db(subject,
                           os.path.join(os.path.dirname(subject),
                                        config.OUTPUT),
                           config)

        # Blast'em
        outputs.append(blastp(query, db, config))

    return outputs


def listify_blast_line(match_line, casts=[]):
    """Return a list with the blast output."""

    new_line = match_line.split("\t")
    for cast in casts:
        # Cast the columns in "casts" to types for further comparisons
        new_line[cast[0]] = getattr(
            __builtin__, cast[1])(new_line[cast[0]])

    return new_line


def get_best_from_blast_output(blast_output):
    """Return only the best matches for each group of matches."""

    blast_list = [listify_blast_line(y, casts=[(IDENTITY, "float")])
                  for y in blast_output.rstrip().split("\n")]

    return simplify_blast_output(blast_list=blast_list, group=[])


def get_best_from_group(blast_list, position):
    """Return the best item from a list."""

    # Sort the matches by position and get the best of them
    best_match = sorted(
        blast_list, key=itemgetter(position), reverse=True)[0]

    return [str(x) for x in best_match]


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
