### Energy Evaluation of Single Mutant Homology Models:
![photo](Media/Picture1.png)
### Installation and Dependencies:
1. Clone repo:
```sh
   git clone https://github.com/eved1018/InterfaceExtention
```
2. Download Modeller: https://salilab.org/modeller/download_installation.html
    * for Conda enviroment:
```sh 
        conda config --add channels salilab
        conda install modeller
```
3. Download python dependencies (pyhull, scipy, numpy):
```sh
    pip install -r requirements.txt 
```
4. Install and compile EvoEf1
```sh
    cd src/
    git clone https://github.com/tommyhuangthu/EvoEF.git
    cd EvoEF
    g++ -O3 -ffast-math -o EvoEF src/*.cpp
```

5. download foldx (https://foldxsuite.crg.eu/) and move the executable and rotobase.txt to src/foldx folder.

### Configuration File:
The config.txt is used to guide the mutagensis and is organzied in trhee columns:
1) `3 letter amino acid` 
2) `residue number`
3) `comma seprated list of amino acids to mutate to or * for all`
if config.txt is not empty you will be prompted to select from the interface positions. 

### Command line Arguments:
* `-pdb`: RCSB PDB id, if not provided you will be prompted to select one. If it is is in the input/ folder it will be used. Otherwise it will be downloaded from the RCSB.
* `-qc`: Query chain to mutate.
* `-ic`: partner chain.