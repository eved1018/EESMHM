import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from .intercaatInterface import *
from .CLI import cli
from .mutant_model import Mutate
from .energy_eval import EnergyCalculation, WTEnergyCalc
from .input_parser import configManager
from .AminoAcidConverter import *
from .heatplot import heatmap


def main(pdb, pdb_file, qc, ic, config,foldx_path, evoef_path):
    # wildtype analysis:
    intercaat_wt_scores : dict = {} # {position: [intercaat_wt_score]}
    intercaat_result = intercaatWT(pdb_file, qc, ic)
    # get positions to mutate
    positions_and_mutations, foldx_path, evoef_path = configManager(config, intercaat_result, foldx_path, evoef_path)

    # get wildtype energy values
    evo_wt_score, foldx_wt_score, wt_gbsa = WTEnergyCalc(pdb_file, qc, ic, evoef_path, foldx_path)

    #perform mutation:
    mutant_folder = f"output/{pdb}/mutants/"
    mutation_jobs: list = []
    jobs : list = []
    results : dict = {}
    for position in positions_and_mutations:
        key = position.replace("_", "")
        interactions = intercaat_result[key][2]
        intercaat_wt_scores[position] = interactions
        mutant_self = key + three_to_one(key[:3])
        results[mutant_self] = [foldx_wt_score, evo_wt_score, 0, interactions, 0] #check this works
        mutants = positions_and_mutations[position]
        for mutantAA in mutants:
            mutation_jobs.append([mutantAA, position, mutant_folder, pdb_file, qc])


    return_vals = MultiThreadedExe(mutation_jobs, Mutate)
    for i in return_vals:
        if len(i) == 1:
            results[i[0]] = ["nan","nan","nan","nan","nan"]
        else:
            position, mutant_name, wt_AA,respos, mutantfile, mutantAA_3letter = i
            jobs.append([position, mutant_name, wt_AA,respos, mutantfile,foldx_wt_score, evo_wt_score,intercaat_wt_scores, mutant_folder, qc, ic, mutantAA_3letter, pdb_file,wt_gbsa, evoef_path, foldx_path])


    return_vals = MultiThreadedExe(jobs, EnergyCalculation)
    for i in return_vals:
        mutant_name, foldx, evoef, intercaat,contacts,md_energy = i
        results[mutant_name] = [foldx, evoef, intercaat,contacts, md_energy]

    df = dataframeWrapper(results)
    df.to_csv(f"output/{pdb}/{pdb}_results.csv")
    heatmap(df,pdb)
    print(df)
    return df

def dataframeWrapper(results):
    df = pd.DataFrame.from_dict(results, orient='index', columns = ["DDG_foldx","DDG_evoef", "Dintercaat_normalized","mutantContacts", "DDG_GBSA" ])
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'mutant'})
    return df


def MultiThreadedExe(jobs, funct):
    with ProcessPoolExecutor() as exe:
        return_vals = exe.map(funct, jobs)
    return return_vals

def EESMHM():
    pdb, pdb_file, qc, ic, config, foldx, evoef = cli()
    df = main(pdb,pdb_file, qc, ic, config, foldx, evoef)
    return
