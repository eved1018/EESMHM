import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from .intercaatInterface import intercaatMutant, intercaatWT
from .CLI import cli
from .mutant_model import mutateModel
from .energy_eval import EvoEF_run, Foldx_run, GBSA
from .input_parser import inputParser
from .AminoAcidConverter import three_to_one, one_to_three, all_aa_one
from .heatplot import heatmap


"""
TODO:

redo all the keys, and amino acid stuff, do not mix three letter and one letter code in main, hand that off to intercaat and modeller interface

add scrwl4

speed up
"""


def main(pdb, pdb_file, qc, ic, config,foldx_path, evoef_path):
    os.makedirs(f"output/{pdb}/mutants/", exist_ok=True)
    os.makedirs(f"tmp/", exist_ok=True)

    #wildtype analysis:
    intercaat_wt_scores : dict = {} # {position: [intercaat_wt_score]}
    intercaat_result = intercaatWT(pdb_file, qc, ic)
    if evoef_path:
        evo_wt_score = EvoEF_run(pdb_file, qc, ic, evoef_path)
    else:
        evo_wt_score = 0
    if foldx_path:
        foldx_wt_score = Foldx_run(pdb_file, qc, ic, foldx_path)
    else:
        foldx_wt_score = 0
    wt_gbsa   = GBSA(pdb_file, qc, ic )
    mutant_folder = f"output/{pdb}/mutants/"

    positions_and_mutations, foldx_path, evoef_path = configManager(config, intercaat_result, foldx_path, evoef_path)

    #perform mutation:
    jobs : list = []
    for position in positions_and_mutations:
        key = position.replace("_", "")
        interactions = intercaat_result[key][2]
        intercaat_wt_scores[position] = interactions

        mutants = positions_and_mutations[position]
        for mutantAA in mutants:
            # mutantAA_3letter = aa_three_letter[mutantAA]

            mutantAA_3letter = one_to_three(mutantAA)
            # create mutant PDB
            wt_AA , respos = position.split("_")
            mutposition = respos + mutantAA
            mutantfile =  mutant_folder + wt_AA + mutposition + ".pdb"
            mutantfile = mutateModel(pdb_file, respos, mutantAA_3letter, qc, mutantfile)
            mutant_name = wt_AA + respos + mutantAA
            jobs.append([position, mutant_name, wt_AA,respos, mutantfile,foldx_wt_score, evo_wt_score,intercaat_wt_scores, mutant_folder, qc, ic, mutantAA_3letter, pdb_file,wt_gbsa, evoef_path, foldx_path])


    # for each mutant run foldx, evoef and intercaat in parrelel
    results : dict = {}
    with ProcessPoolExecutor() as exe:
        return_vals = exe.map(mutantEnergyScorer, jobs)
        for return_val in return_vals:
            mutant_name, foldx, evoef, intercaat,contacts,md_energy = return_val
            results[mutant_name] = [foldx, evoef, intercaat,contacts, md_energy]

    # convert dict to dataframe
    df = pd.DataFrame.from_dict(results, orient='index', columns = ["DDG_foldx","DDG_evoef", "Dintercaat_normalized","mutantContacts", "DDG_GBSA" ])
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'mutant'})
    print(df)
    print(pdb, pdb_file)
    df.to_csv(f"output/{pdb}/{pdb}_results.csv")
    heatmap(df)
    return df

def mutantEnergyScorer(params):
    evoef = 0
    foldx = 0
    position, mutant_name, wt_AA, respos, mutantfile,foldx_wt_score, evo_wt_score,intercaat_wt_scores, mutant_folder, qc, ic,mutant_AA, pdb_file , wt_gbsa, evo_path, foldx_path = params
    # print(mutant_name)
    if evo_path:
        evoef = EvoEF_run(mutantfile, qc, ic, evo_wt_score, evo_path)
    mutant_pdb = mutantfile.split("/")[-1]
    if foldx_path:
        foldx = Foldx_run(mutant_pdb, qc,ic, mutant_folder,foldx_wt_score, foldx_path)
    dintercaat, contacts = intercaatMutant(position, mutant_pdb, qc, ic, mutant_folder, intercaat_wt_scores, wt_AA, mutant_AA )
    md_energy = GBSA(mutantfile, qc, ic, wt_gbsa )
    return mutant_name, foldx, evoef, dintercaat, contacts, md_energy

def configManager(config, intercaat_result, foldx_path, evoef_path):
    positions_and_mutations, foldx_path, evo_path = inputParser(config, foldx_path, evoef_path)
    if len(positions_and_mutations) == 0:
        positions_and_mutations = {}
        for i in list(intercaat_result.keys()):
            wt_aa_three = i[:3]
            i = i[:3] + "_" + i[3:]
            positions_and_mutations[i] = [i for i in all_aa_one() if one_to_three(i) != wt_aa_three]
    return positions_and_mutations, foldx_path, evo_path


def EESMHM():
    pdb, pdb_file, qc, ic, config, foldx, evoef = cli()
    df = main(pdb,pdb_file, qc, ic, config, foldx, evoef)
    return
