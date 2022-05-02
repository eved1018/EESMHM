from intercaat.intercaatWrapper import intercaat
from halo import Halo

@Halo(text='Loading', spinner='dots')
def intercaatWT(pdb_file, qc, ic):
    match,  interactions = intercaat(pdb_file, qc, ic)
    return interactions

def intercaatMutant(position, pdb_file, qc, ic, fp, intercaat_wt_scores, wt_AA, mutant_AA ):

    avergae_interactions = {
                "ALA":5.054054054054054,
                "APR":6.0,
                "ARG":10.886363636363637,
                "ASN":7.120481927710843,
                "ASP":7.112903225806452,
                "CYS":5.45,
                "GLN":7.078651685393258,
                "GLU":9.020689655172413,
                "GLY":6.009174311926605,
                "HIS":9.796610169491526,
                "ILE":7.917355371900826,
                "LEU":8.095238095238095,
                "LYS":8.581560283687944,
                "MET":8.830508474576272,
                "PHE":13.96124031007752,
                "PRO":6.953125,
                "SER":7.407692307692308,
                "THR":6.0,
                "TRP": 16.45,
                "TYR":15.330645161290322,
                "VAL":6.371621621621622,}
    key = mutant_AA + position.split("_")[1]
    match, results = intercaat(pdb_file,qc,ic,mi=0,fp= fp)
    if key in results:
        interactions = results[key][2]
    else:
        interactions = 0
    interactions_normalized = interactions/avergae_interactions[mutant_AA]
    wt_interactions = intercaat_wt_scores[position] / avergae_interactions[wt_AA]
    dinteractions  = wt_interactions - interactions_normalized
    return round(dinteractions, 3), interactions




