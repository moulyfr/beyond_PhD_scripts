These 2 scripts are used to run sc velocity analysis with scvelo.

If working with Seurat objects, then the first script needs to be run to extract raw, normalized counts and metadata.

Then, scVelo is performed with the python script. There are 3 different models that can be run, briefly:
1. Stochastic model: probabilistic
- considers noise and variability
- use when data not well defined
2. Dynamic model: deterministic
- fits continuous time dynamics
- use when dynamics well defined or time series data avail
3. Kinetic model : integrated
- accounts for both spliced+ unspoiled
- use when understanding of RNA prod/degradation needed
