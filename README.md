# EMTanalysis
EMT transcriptional and metabolomic analysis

This repository contains the code and input files for the analysis of Non-small cell lung cancer cell lines data. 

Transcriptional and metabolomic data was retrieved from  [Sun Y., *et. al.* 2014](https://doi.org/10.1186/2049-3002-2-20).

We used differential expression analysis for transcriptional data. Using frma and limma R packages.

Besides, we used an algorithm called [Dycone](https://github.com/cdiener/dycone) for metabolomic analysis. 
Briefly, the algorithm assumes that all the metabolic reactions in reconstruction obey the law-mass action and the system has reached a steady-state condition.
Under these constraints, the feasible space of metabolic phenotypes before and after the EMT can be obtained.

Combining these results with microarray analysis we were able to identify specific reactions that change in three different NSCLC cell lines with different genetic background.

