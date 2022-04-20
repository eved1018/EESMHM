import subprocess

def EvoEF_run(pdb, qc, ic, wt_evo_score = 0.0):
    cmd = f"src/EvoEF/EvoEF --command=ComputeBinding --split={qc},{ic} --pdb={pdb}|tail -n2"
    p = subprocess.run(cmd, shell = True,check = True, stdout=subprocess.PIPE, universal_newlines=True)
    try:
        score = p.stdout.split()[2]
        score =float(score) - wt_evo_score
        score = round(score, 3)
    except:
        score = "nan"
    return score
        
def Foldx_run(mutant, qc,ic, mutant_folder, wt_foldx_score = 0.0):
    out = subprocess.run(f"./src/Foldx/foldx --command=AnalyseComplex --analyseComplexChains={qc},{ic} --pdb={mutant} --pdb-dir={mutant_folder} --output-dir=tmp/ |tail -n 7|head -n 1 ", shell = True,stdout=subprocess.PIPE, universal_newlines = True)
    out = out.stdout
    out = str(out)
    mutant_dg = float(out.split()[2])
    score = mutant_dg - wt_foldx_score
    return round(score, 3) 
