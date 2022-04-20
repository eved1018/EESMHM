import os
from concurrent.futures import ProcessPoolExecutor
from unittest import result
import pandas as pd
import inquirer

from src.intercaat import intercaatMutant, intercaatWT
from src.CLI import cli
from src.mutant_model import mutateModel
from src.energy_eval import EvoEF_run, Foldx_run
from src.input_parser import inputParser



def main(pdb, pdb_file, qc, ic, config):
    aa_three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 
                    'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    
                    'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',  
                    'G':'GLY', 'P':'PRO', 'C':'CYS'}

    os.makedirs(f"output/{pdb}/mutants/", exist_ok=True)
    os.makedirs(f"tmp/", exist_ok=True)

    #wildtype analysis:
    intercaat_wt_scores : dict = {} # {position: [intercaat_wt_score]}
    intercaat_result = intercaatWT(pdb_file, qc, ic)
    evo_wt_score = EvoEF_run("input/"+pdb_file, qc, ic)
    foldx_wt_score = Foldx_run(pdb_file, qc, ic, "input/")
    mutant_folder = f"output/{pdb}/mutants/"

    positions_and_mutations = configManager(config, intercaat_result)

    #perform mutation:
    jobs : list = []
    for position in positions_and_mutations:

        key = position.replace("_", "")
        interactions = intercaat_result[key][2]
        intercaat_wt_scores[position] = interactions

        mutants = positions_and_mutations[position]
        for mutantAA in mutants: 
            mutantAA_3letter = aa_three_letter[mutantAA]
            # create mutant PDB
            wt_AA , respos = position.split("_")
            mutposition = respos + mutantAA 
            mutantfile =  mutant_folder + wt_AA + mutposition + ".pdb"
            mutantfile = mutateModel(pdb_file, respos, mutantAA_3letter, qc, mutantfile, "input/")
            mutant_name = wt_AA + respos + mutantAA
            jobs.append([position, mutant_name, wt_AA,respos, mutantfile,foldx_wt_score, evo_wt_score,intercaat_wt_scores, mutant_folder, qc, ic, mutantAA_3letter, pdb_file])

    
    # for each mutant run foldx, evoef and intercaat in parrelel
    results : dict = {}
    with ProcessPoolExecutor() as exe:
        return_vals = exe.map(mutantEnergyScorer, jobs)
        for return_val in return_vals:
            mutant_name, foldx, evoef, intercaat = return_val
            results[mutant_name] = [foldx, evoef, intercaat]

    # convert dict to dataframe 
    df = pd.DataFrame.from_dict(results, orient='index', columns = ["DDG_foldx","DDG_evoef", "Dintercaat_normalized" ])
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'mutant'})
    print(df)
    df.to_csv(f"output/{pdb}/{pdb}_results.csv")
    return 

def mutantEnergyScorer(params):
    position, mutant_name, wt_AA, respos, mutantfile,foldx_wt_score, evo_wt_score,intercaat_wt_scores,  mutant_folder, qc, ic,mutant_AA, pdb_file  = params
    evoef = EvoEF_run(mutantfile, qc, ic, evo_wt_score)
    mutant_pdb = mutantfile.split("/")[-1]
    foldx = Foldx_run(mutant_pdb, qc,ic, mutant_folder,foldx_wt_score)
    dintercaat = intercaatMutant(position, mutant_pdb, qc, ic, mutant_folder, intercaat_wt_scores, wt_AA, mutant_AA )
    return mutant_name, foldx, evoef, dintercaat

def configManager(config, intercaat_result):
    amino_acids = ['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
    positions_and_mutations = inputParser(config)
    if len(positions_and_mutations) == 0:
        questions = [
            inquirer.Checkbox('mutants',
                    message="Select which positions to mutate (using left arrow)?",
                    choices=list(intercaat_result.keys()),
                    ),
        ]
        answers = inquirer.prompt(questions)
        answers = list(answers.values())[0]
        positions_and_mutations = {}
        for i in answers:
            i = i[:3] + "_" + i[3:]
            positions_and_mutations[i] = amino_acids
    return positions_and_mutations


if __name__ == "__main__":
    pdb, pdb_file, qc, ic, config = cli()
    main(pdb, pdb_file, qc, ic, config)