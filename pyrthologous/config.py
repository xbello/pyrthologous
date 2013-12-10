"""Config the user constants."""
import os

MUSCLE = os.path.expanduser("~/bin/muscle")

BLAST_PATH = "/usr/bin/"
BLASTP = os.path.join(BLAST_PATH, "blastp")
BLAST_DB_MAKER = "makeblastdb"

SUFFIX = "NORM"
SEP = ";"
BASE_PATH = os.path.join(os.getcwd())
LEAF_PATHS = [
    "Babesia",
    "Giardia",
    "Theileria",
    "Toxo"
]


# This ORGANISMS dict relates a "real" genome filename with a shorting of it.
# I.e. for a fasta called "Babesia_bovis.fas" we have a short of "Bv" in the
# form of key:value dict like {"Bv": "Babesia_bovis"}
ORGANISMS = {
    "Bv": "Babesia_bovis",
    "Bm": "Babesia_microti",
    "GA": "Giardia_A",
    "GB": "Giardia_B",
    "GE": "Giardia_E",
    "Tp": "Theileria_parva",
    "Ta": "Theileria_annulata",
    "TG": "Toxo_GT",
    "TM": "Toxo_ME",
    "TV": "Toxo_VE"}

OUTPUT = "outputs"
