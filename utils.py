import itertools
import os
from config import BASE_PATH

def make_reciprocals(src_dir):
    """String -> List

    Where String is the name of the leaf which contains all genomes and
    List is the list of permutations to be generated
    """
    genome_src = os.path.join(BASE_PATH, src_dir)
    reciprocals = [x for x in os.listdir(genome_src) if x.endswith("fas")]
    permuts = [x for x in itertools.combinations(reciprocals, 2)]

    return permuts

