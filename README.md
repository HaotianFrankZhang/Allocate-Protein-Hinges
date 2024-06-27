# Allocate-Protein-Hinges
Codes and method details of the article "Hinge sites of proteins are alternative target sites for drug binding".

#  Prepare the packages and jupyter notebooks

Some packages are needed. 

```terminal
pip install numpy
pip install prody
pip install matplotlib
pip install scipy
pip install notebook
```
If you use annaconda vitual environment, you can use "environment.yml" to buid a environment to run the code.
```terminal
conda env create -f environment.yml -n hinges
conda activate hinges
```

# Program details and methods (See Work_flow.ipynb)
### implement pacakges
To import all functions, you need to download the "GNMhinges.py" in the main folder. Please make sure it exists in the same directory when you running the program
```python
from numpy import *
from prody import *
from matplotlib.pyplot import *
from scipy.stats import hypergeom
import os
from GNMhinges import *
```

### Build ensemble by query Dali
Suppose you are interested in B2AR (PDB: 3d4s chain A)
Similar ensemble, length difference > 0.8; RMSDs less than 2; Z score greater that 10
```python
eigVals, averageEigVects, ids, gnms = getModesSimilar('3d4s', 'A', length=0.8, rmsd=2, Z=10)
```
Diverse ensemble, length difference <= 0.95; RMSDs greater than 1; Z score greater that 10
```python
eigVals, averageEigVects, ids, gnms = getModesDiverse('3d4s', 'A', length=0.95, rmsd=1, Z=10)
```
Ensemble for dimer, you need to query Dali twice for each chains, and build your ensemble based on the result ids. Take COX2 as an example
```python
currPDB = '4m11'
eachChain = 'A'
averageEigVals_A, averageEigVects_A, ids_A, gnms_A = getModesDiverse(currPDB, eachChain, length=0.95, rmsd=1, Z=10)
eachChain = 'B'
averageEigVals_B, averageEigVects_B, ids_B, gnms_B = getModesDiverse(currPDB, eachChain, length=0.95, rmsd=1, Z=10)
```
### Optional: set up the reference structure
If you are only interested in part of the structure. For B2AR, we are only interested in 7-TM domain. The rest part in the crystal structure is just used to stablize the whole structure
```python
GPCR = parsePDB('3d4s', subset='calpha') 
calphas = GPCR.select('ca resnum 1:343') # your reference structure, not whole protein but only 7-transmemebrane domain
```
Note: The default reference structure is first element in your ids list.

### Parse your ensembles
Parse the PDB files for the specified IDs. Here we focus on the alpha carbons (`ca`) to simplify the model.
```python
ags = parsePDB(ids_similar, subset='ca')
```
Build ensembles where the reference structure by default is the first item. You can specify a different reference structure by setting the ref parameter.
```python
dali_ens = buildPDBEnsemble(ags, ref=calphas)
```
Perform calculations on the ensemble to get normal mode analysis using the Gaussian Network Model (GNM).
```python
gnms = calcEnsembleENMs(dali_ens, model='GNM', trim='reduce', n_modes=None)
```
Extract eigenvalues and eigenvectors from the GNM analysis.
```python
eigVals = gnms.getEigvals()
averageEigVals = gnms.getEigvals()[0]
eigVects = gnms.getEigvecs()
averageEigVecs = mean(eigVects, axis=0)
```
### Num. of modes are considered to allocate hinges
It is based on the EigenValues. You can change the contribution threshold by changing the second parameter
```python
currNumModes = getModesGivenThreshold(averageEigVals, 0.33)
```

### Allocate hinges from GNM modes.
You need average Eigenvectors. You can change the number of modes by changing the second parameter. You can decrease the third parameter to adjust the band to include more residues as hinges.
```python
Hinges_3modes = getHinges2(averageEigVecs, 3, 15) # using 3 modes
```
Optional:The function trimEnds is designed to refine hinge predictions by removing residues at the protein ends that are unlikely to function as hinges. This function takes three parameters:
1. Raw Hinges Indices (Hinges_3modes): This parameter contains the indices of potential hinge points identified in the protein structure.
2. Protein Chain Boundaries: This list specifies the start and end indices of each protein chain or subunit within the dimer or complex. For example, if the reference protein is a dimer where the first monomer's residues range from 0 to 551 and the second monomer's residues range from 552 to 956, the parameter should be provided as [[0, 551], [552, 956]]. Each subunit or chain in the protein is represented by a list, where the first element is the starting residue index, and the second element is the ending residue index.
3. Exclusion Zone (trimEnds): This integer specifies how many residues at the beginning and end of each chain should be excluded from consideration as hinge points, as they are less likely to act as hinges due to their position.
For example, when calling the function with:
```python
Hinges_3modes_final = trimEnds(Hinges_3modes, [[0, 278]], 20)
```
### Define binding sites.
Define the binding sites. Please make sure the binding sites are indeces of the residues, which are always starting from 0. In B2AR, we found 37 binding sites:
```python
bindings = [27,38,41,42,45,48,49,52,53,76,77,78,80,81,82,83,85,86,119,122,123,126,134,161,163,167,168,171,172,175,222,225,226,229,244,248,252]
```
Make sure all indeces are integers
```python
binding = list(set([int(x) for x in bindings]))
```

### Calculate p values.
N: total number of residues in your protein <br>
```python
protein_length = 279
```
M: # of drug binding residues <br>
```python
M = len(binding)
```
n: # of hinge residues <br>
```python
n = len(Hinges_3modes_final)
```
k: overlaps 
```python
overlaps_2 = len(binding) + len(Hinges_3modes_final) - len(set(binding + Hinges_3modes_final))
```
p values
```python
HyperScore_2 = ORA(len(binding), protein_length, len(Hinges_3modes_final), overlaps_2) # p_value = ORA(M, N, n, k)
```
Print the final results
```python
print ('# of binding sites is', len(binding))
print ('# of hinge sites for threshold 0.33, overlap, hyper score', len(Hinges_3modes_final), overlaps_2, HyperScore_2)
```

### Other functions and more details
You can find other functions in "Allocate-Protein-HInges.ipynb". For the details of each case (including results for both similar and diverse ensemble), you can find jupyter notebook files in 20 folders named by the protein acronyms.

To plot the final graphs, you need the mode you want to plot, average Eigenvectors, gnms, binding sites. If you want to mark other sites you are interested in, you can file the index in the last parameter
```python
modes = 0 # plot first mode
plotSingleGraph(modes, averageEigVecs, gnms, binding, [])
```
In the program, many functions have been implemented in prody pacakges. To check the details of those functions, please visit the prody website to explore more: http://www.bahargroup.org/prody/

# Raw data
Within each of the 20 folders, there is a subfolder named "Data" that contains several key files: the ensemble ID list, Dali search results, and information on binding sites. It is important to note that the results from Dali searches may vary over time since the Protein Data Bank (PDB) database is continually updated. To ensure consistency in your data, it is essential to maintain records of the ensemble lists and the corresponding chain IDs.

Additionally, located in the main folder, you will find a file called "drug_letters_screen.txt". This file comprises a list of 1,018 approved drug IDs, providing a comprehensive resource for further analysis.

