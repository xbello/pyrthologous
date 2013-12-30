"""Config the user constants."""
import os

MUSCLE = os.path.expanduser("~/bin/muscle")

BLAST_PATH = "/usr/bin/"
BLAST_PATH = "C:\\Program files\\NCBI\\blast-2.2.28+\\bin"
BLASTP = os.path.join(BLAST_PATH, "blastp.exe")
BLAST_DB_MAKER = os.path.join(BLAST_PATH, "makeblastdb.exe")

SUFFIX = "NORM"
SEP = ";"
BASE_PATH = os.path.join(os.getcwd())
# This is where the genomes are relative to BASE_PATH:
GENOMES = "mock_genome"

OUTPUT = "outputs"
