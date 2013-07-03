#!/usr/bin/env python
import os

base_path = "/opt/space/data/Parasitos/secuencias"
src_dirs = ["Babesia", "Giardia", "Theileria", "Toxoplasma"]

for src_dir in src_dirs:
    genome_src = os.path.join(base_path, src_dir)
    for fasta in os.listdir(genome_src):
        if fasta.endswith("fasta"):
            os.system("./rehead_fasta.py {0} {1} -t".format(
                os.path.join(genome_src, fasta),
                os.path.join(genome_src, fasta.split(".")[0] + ".txt")))
