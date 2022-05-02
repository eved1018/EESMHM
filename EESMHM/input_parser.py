import os
from .AminoAcidConverter import *

def inputParser(input_file, foldx_path, evoef_path):
    positions_and_mutations = {}
    if os.path.exists(input_file):
        with open(input_file, "r") as infile:
            for line in infile.readlines():
                if line.startswith("#foldx") and len(line.split()) > 1:
                    foldx_path = line.split()[1]
                if line.startswith("#evoef") and len(line.split()) > 1:
                    evoef_path = line.split()[1]
                wt_aa, position, mutations = line.split()
                if len(mutations) == 1:
                    if mutations == "*":
                        mutations = [i for i in all_aa_one() if i != wt_aa ]
                    else:
                        mutations = list(mutations)
                else:
                    mutations = mutations.split(",")
                positions_and_mutations[one_to_three(wt_aa)+ "_" +position] = mutations
    return positions_and_mutations, foldx_path, evoef_path
