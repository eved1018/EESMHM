### Energy Evaluation of Single Mutant Homology Models:

### Installation and Dependencies:
1. Install:
```sh
   pip install EESMHM
```
2. Download Modeller: https://salilab.org/modeller/download_installation.html
    * for Conda enviroment:
```sh 
        conda config --add channels salilab
        conda install modeller
```
3. Download EvoEF (https://github.com/tommyhuangthu/EvoEF).

4. Download Foldx (https://foldxsuite.crg.eu/).


### Configuration File:
The config.txt is used to guide the mutagensis and is organzied in two parts

part 1: paths
```sh
	#foldx {foldx path}
	#evoef {evoef path}
```

part 2: mutation config:
`1 letter amino acid code` `residue number` `comma seprated list of amino acids to mutate to or * for all`

If left empty all interface positions will be mutated.  

### Command line Arguments:
* `-pdb`: RCSB PDB id, if not provided you will be prompted to select one. If it is is in the current working directory it will be used. Otherwise it will be downloaded from the RCSB.
* `-qc`: Query chain to mutate.
* `-ic`: Partner chain.
* `-config`: config file path
* `-foldx`: foldx path
* `-evoef`: Evoef path
