# These various scripts conduct the following with processed seurat objects:

1) cell-cell communication analysis with CellPhoneDB (cpdb) pkg: cpdb-related input prep, running the analysis, filtering results, and lastly conducting a meta-analysis
   
2) cell-cell communication analysis with CellChat pkg to run initial analysis

3) gene ontology (GO) analysis

4) disease vs healthy cell proportion analysis

5) disease vs healthy differential gene expression (DEG) analysis (with MAST pkg)

6) sc velocity analysis (using scVelo pkg)
-  in general terms, is used to extract genes involved during dynamic processess (of cell development, differentiation, etc) by examining unspliced vs spliced reads
-  as both unspliced and spliced reads are seq'd (and such info is avaialble in the bam files), these are used to infer time, since right after a gene turns on, the cell will have unspliced transcripts
-  high velocity of a particular cell = gene was just turned on 

7) general graphing scripts
