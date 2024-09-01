## These 2 scripts are used to run sc velocity analysis with scvelo.

1) First, the bam files are converted to loom files using Velocyto (shell script).

2) Then, scVelo is performed with the python script (after the loom is converted to anndata). 

## There are 3 different models that can be run, briefly:

### A. Stochastic model: probabilistic
- considers noise and variability
- use when data not well defined

### B. Dynamic model: deterministic
- fits continuous time dynamics
- use when dynamics well defined or time series data avail

### C. Kinetic model : integrated
- accounts for both spliced+ unspoiled
- use when understanding of RNA prod/degradation needed
