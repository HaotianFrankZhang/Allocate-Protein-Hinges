# Allocate-Protein-Hinges
Codes and method details of the article "Hinge sites of proteins are alternative target sites for drug binding".

#  Prepare the packages and jupyter notebooks

Some packages are needed. The details are in "environment.yml"

```terminal
pip install numpy
pip install prody
pip install matplotlib
pip install scipy
pip install notebook
```

# Code explanation

### Build ensemble by query Dali
Suppose you are interested in B2AR (PDB: 3d4s chain A)
```python
# Similar ensemble, length difference <= 0.8; RMSDs less than 2; Z score greater that 10
eigVals, averageEigVects, ids, gnms = getModesSimilar('3d4s', 'A', length=0.8, rmsd=2, Z=10)
# Diverse ensemble, length difference <= 0.95; RMSDs greater than 1; Z score greater that 10
eigVals, averageEigVects, ids, gnms = getModesDiverse('3d4s', 'A', length=0.95, rmsd=1, Z=10)
```

### Parse your ensembles
Parse the PDB files for the specified IDs. Here we focus on the alpha carbons (`ca`) to simplify the model.
```python
ags = parsePDB(ids_similar, subset='ca')
```
Build ensembles where the reference structure by default is the first item. You can specify a different reference structure by setting the ref parameter.
```python
dali_ens = buildPDBEnsemble(ags, ref=calphas)
```
### Calculate Ensemble ENMs
Perform calculations on the ensemble to get normal mode analysis using the Gaussian Network Model (GNM).
```python
gnms = calcEnsembleENMs(dali_ens, model='GNM', trim='reduce', n_modes=None)
```
### Retrieve Eigenvalues and Eigenvectors
Extract eigenvalues and eigenvectors from the GNM analysis. Additionally, compute averages as needed.
```python
eigVals = gnms.getEigvals()
averageEigVals = gnms.getEigvals()[0]
eigVects = gnms.getEigvecs()
averageEigVecs = mean(eigVects, axis=0)
```
### How many modes are needed to reach 33% contribution based on the Eigenvalues.
You can change the contribution threshold by changing the second parameter
```python
currNumModes = getModesGivenThreshold(averageEigVals, 0.33)
```
### Allocate hinges from GNM modes.
You need average Eigenvectors. You can change the number of modes by changing the second parameter. You can decrease the third parameter to adjust the band to include more residues as hinges.
```python
Hinge_index = getHinges2(averageEigVecs, 3, 15)
```

### Calculate p values.
N: total number of residues in your protein <br>
M: # of drug binding residues <br>
n: # of hinge residues <br>
k: overlaps 
```python
p_value = ORA(M, N, n, k)
```

### Plot graphs
You need the mode you want to plot, average Eigenvectors, gnms, binding sites. If you want to mark other sites you are interested in, you can file the index in the last parameter
```python
modes = 0 # plot first mode
plotSingleGraph(modes, averageEigVecs, gnms, binding, [])
```

