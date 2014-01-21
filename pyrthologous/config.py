"""Config the user constants."""
import os

MUSCLE = os.path.expanduser("~/bin/muscle")

BLAST_PATH = "/usr/bin/"
BLASTP = os.path.join(BLAST_PATH, "blastp")
BLAST_DB_MAKER = os.path.join(BLAST_PATH, "makeblastdb")

SUFFIX = "NORM"
SEP = ";"
BASE_PATH = os.path.join(os.getcwd())
# This is where the genomes are relative to BASE_PATH:
GENOMES = "Parasites"

OUTPUT = "outputs"
