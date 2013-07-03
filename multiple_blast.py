#!/usr/bin/env python
import itertools
import os

base_path = "/opt/space/data/Parasitos/secuencias"
src_dirs = ["Babesia", "Giardia", "Theileria", "Toxoplasma"]

def make_output_name(file_1, file_2):
    """String, String -> String

    Takes two filenames (standarized) and return the output filename
    """
    #This takes a name like "Babesia_microti_NORM.fas" and get a "B_m"
    query = "".join([x[0] for x in file_1.split("_")[:2]])
    subject = "".join([x[0] for x in file_2.split("_")[:2]])

    return "{0}_vs_{1}.blast".format(query, subject)

for src_dir in src_dirs:
    genome_src = os.path.join(base_path, src_dir)
    reciprocals = [x for x in os.listdir(genome_src) if x.endswith("fas")]
    permuts = [x for x in itertools.combinations(reciprocals, 2)]

    for permut in permuts:
        for x, y in [(0, 1), (1, 0)]:
            output = make_output_name(permut[x], permut[y])
            os.system("python reciprocal_blast.py {0} {1} > {2}".format(
                os.path.join(genome_src, permut[x]),
                os.path.join(genome_src, permut[y]),
                os.path.join(os.getcwd(), "outputs", output)))
