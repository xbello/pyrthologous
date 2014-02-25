"""Config the user constants."""
import logging
import os

MUSCLE = os.path.expanduser("~/bin/muscle")

BLAST_PATH = "/usr/bin/"
BLASTP = os.path.join(BLAST_PATH, "blastp")
BLAST_DB_MAKER = os.path.join(BLAST_PATH, "makeblastdb")

BASE_PATH = os.path.join(os.getcwd())
LOGFILE=os.path.join(BASE_PATH, "logfile.log")
LOGGER = logging.basicConfig(filename=LOGFILE,
                             filemode='w',
                             format='%(asctime)s %(message)s',
                             level=logging.INFO)

# This is where the genomes are relative to BASE_PATH:
GENOMES = "Parasites"

OUTPUT = "Parasites_output"

COMPARE = (
    ("Ehellem.fasta", "Eromaleae.fasta"),
    ("Ecuniculi.fasta", "Eintestinalis.fasta"),
    ("Edispar.fasta", "Emoshkovskii.fasta"),
    ("Torientalis.fasta", "Tequi.fasta"),
    ("Ehistolytica.fasta", "Enuttalli.fasta"),
    ("Edispar.fasta", "Ehistolytica.fasta"),
    ("Lmajor.fasta", "Lmexicana.fasta"),
    ("Einvadens.fasta", "Emoshkovskii.fasta"),
    ("Ltarentolae.fasta", "Lbrazilensis.fasta"),
    ("Linfantum.fasta", "Ldonovani.fasta"),
)
