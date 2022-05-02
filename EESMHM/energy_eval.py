import subprocess
from modeller import *
from modeller import gbsa
from modeller.scripts import complete_pdb


def EvoEF_run(pdb, qc, ic, wt_evo_score = 0.0, evoef_path = ""):
    cmd = f"{evoef_path}/EvoEF --command=ComputeBinding --split={qc},{ic} --pdb={pdb}|tail -n2"
    p = subprocess.run(cmd, shell = True,check = True, stdout=subprocess.PIPE, universal_newlines=True)
    try:
        score = p.stdout.split()[2]
        score =float(score) - wt_evo_score
        score = round(score, 3)
    except:
        score = "nan"
    return score

def Foldx_run(mutant, qc,ic, mutant_folder="./", wt_foldx_score = 0.0, foldx_path = ""):
    out = subprocess.run(f"{foldx_path}/foldx --command=AnalyseComplex --analyseComplexChains={qc},{ic} --pdb={mutant} --pdb-dir={mutant_folder} --output-dir=tmp/ |tail -n 7|head -n 1 ", shell = True,stdout=subprocess.PIPE, universal_newlines = True)
    try:
        out = out.stdout
        out = str(out)
        mutant_dg = float(out.split()[2])
        score = mutant_dg - wt_foldx_score
        score = round(score, 3)
    except:
        score = "nan"
    return score

def GBSA(pdb,qc: str,ic: str ,wt_gbsa = 0.0):
    log.none()
    env = Environ()
    env.io.atom_files_directory = ['../atom_files']
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Calculate just the GB/SA score; turn off soft-sphere
    env.edat.dynamic_sphere = False
    env.edat.energy_terms.append(gbsa.Scorer())
    # GB/SA falls off slowly with distance, so a larger cutoff than the
    # default (4.0) is recommended
    env.edat.contact_shell = 8.0
    mdl = complete_pdb(env, pdb)

    # Select all atoms
    # atmsel = Selection(mdl)
    atmsel = Selection(mdl)
    atmsel.add(mdl.chains[qc])
    atmsel.add(mdl.chains[ic])
    # Calculate the energy
    energy = atmsel.energy()[0]
    ddg = energy - wt_gbsa
    return ddg
