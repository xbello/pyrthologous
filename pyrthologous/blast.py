"""Wrapper for blast program."""

import __builtin__
import os
import subprocess
from operator import itemgetter
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

    stdout = get_best_from_blast_output(stdout)

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
    """Blast the first of query_subject against the second and viceversa.

    query_subject is a set as (query, subject)

    """

    for query, subject in query_subject, query_subject[::-1]:
        # Make the database
        db = make_blast_db(subject,
                           os.path.join(os.path.dirname(subject), OUTPUT))
        # Blast'em
        stdout, stderr = blastp(query, db)
        # TODO: delete the database
        # Yield and repeat
        yield stdout, stderr


def listify_blast_output(blast_output, casts=[]):
    """Return a list wit the blast output."""

    return_list = []

    for match_line in blast_output.rstrip().split("\n"):
        new_line = match_line.split("\t")
        for cast in casts:
            new_line[cast[0]] = getattr(
                __builtin__, cast[1])(new_line[cast[0]])
        return_list.append(new_line)

    return return_list


def get_best_from_blast_output(blast_output):
    """Return only the best matches for each group of matches."""

    # Identify the columns
    # TODO: These seems to be globals here
    IDENTITY = 2

    blast_list = listify_blast_output(
        blast_output, casts=[(IDENTITY, "float")])

    best_matches = simplify_blast_output(blast_list=blast_list,
                                         group=[],
                                         best_matches=[])

    return "\n".join(["\t".join([str(x) for x in y]) for y in best_matches])


def get_best_from_group(blast_list, position):
    """Return the best item from a list."""

    # Sort the matches by position and get the best of them
    best_match = sorted(
        blast_list, key=itemgetter(position), reverse=True)[0]

    return best_match


def simplify_blast_output(blast_list=[],
                          group=[],
                          best_matches=[]):
    """Return only the best match for each group of matches."""

    # Identify the columns
    QUERY = 0
    IDENTITY = 2

    if not blast_list:
        # Process the last group
        best_matches.append(get_best_from_group(group, IDENTITY))
    else:
        this_line = blast_list.pop(0)
        if not group:
            # Start a new group
            group.append(this_line)
        else:
            # A group has been initialized
            if this_line[QUERY] == group[0][QUERY]:
                # The line is in the same group, add and continue
                group.append(this_line)
            else:
                # The line is the first of a new group
                # Get the best line of the group
                best_matches.append(get_best_from_group(group, IDENTITY))
                # Start a new group
                group = [this_line]

        simplify_blast_output(blast_list=blast_list,
                              group=group,
                              best_matches=best_matches)

    return best_matches
