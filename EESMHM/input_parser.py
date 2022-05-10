import os
from .AminoAcidConverter import *

def configManager(config, intercaat_result, foldx_path, evoef_path):
    positions_and_mutations, foldx_path, evo_path = inputParser(config, foldx_path, evoef_path)
    if len(positions_and_mutations) == 0:
        positions_and_mutations = {}
        for i in list(intercaat_result.keys()):
            wt_aa_three = i[:3]
            i = i[:3] + "_" + i[3:]
            positions_and_mutations[i] = [i for i in all_aa_one() if one_to_three(i) != wt_aa_three]
    return positions_and_mutations, foldx_path, evo_path


def inputParser(input_file, foldx_path, evoef_path):
    positions_and_mutations = {}
    if os.path.exists(input_file):
        print(f"reading configs from {input_file}:")
        with open(input_file, "r") as infile:
            for line in infile.readlines():
                if line.strip():
                    if line.startswith("#foldx") and len(line.split()) > 1:
                        foldx_path = line.split()[1]
                    elif line.startswith("#evoef") and len(line.split()) > 1:
                        evoef_path = line.split()[1]
                    else:
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
