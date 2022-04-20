
def inputParser(input_file):
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
            'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
            'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
            'GLY':'G', 'PRO':'P', 'CYS':'C'}
    amino_acids = ['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
    
    positions_and_mutations = {}
    with open(input_file, "r") as infile:
        for line in infile.readlines():
            wt_aa, position, mutations = line.split()
            if len(mutations) == 1:
                if mutations == "*":
                    wt_aa_3  = one_letter[wt_aa]
                    mutations = [i for i in amino_acids if i != wt_aa_3 ]
                else:
                    mutations = list(mutations)
            else:
                mutations = mutations.split(",")
            positions_and_mutations[wt_aa+ "_" +position] = mutations
    return positions_and_mutations