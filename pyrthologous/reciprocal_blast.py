import os
import subprocess
from config import SUFFIX, BLASTP, BLAST_DB_MAKER


def blastp(query, subject):
    """Make a blastp between query and subject."""

    os.system("{0} -query {1} -db {2} -outfmt 6 -max_target_seqs 1".format(
        BLASTP, query, subject))

    return True

def make_blast_db(tgt):
    '''String -> Bool

    Format a FASTA file as suitable blastp db
    src is the absolute path to the fasta file.'''

    os.system("{0} -in {1} -dbtype prot -logfile DBBLAST.log".format(
        BLAST_DB_MAKER, tgt))

    return True

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("query",
        help="Name of the fasta query file.")
    parser.add_argument("subject",
        help="Name of the fasta subject file")

    args = parser.parse_args()

    file_subject = os.path.join(os.path.abspath(os.getcwd()), args.subject)
    file_query = os.path.join(os.path.abspath(os.getcwd()), args.query)

    make_blast_db(tgt = file_subject)
    blastp(query = file_query, subject = file_subject)
