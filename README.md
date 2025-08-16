# Breast_Cancer_Survival_Analysis_Final_Project
Comparative Analysis to evaluate whether multi-omics integration approaches  such as  MOFA2, DIABLO improve breast cancer survival prediction compared to single-omics  models.

The aim of this study is to compare and evaluate survival prediction performance across single-omics (RNA,CNV) and multi-omics approaches such as MOFA2 and DIABLO using the TCGA-BRCA dataset via UCSC Xena Browser.

## Structure of the project
- 'Breast_Cancer_Survival_Analysis.R' : Complete R script containing data processing, model training, survival analysis and pathway enrichment
- 'results': contains plots (e.g ,Kaplan-Meier plots, pathway enrichment plot, data exploration plots etc



## Dataset Source
The dataset used for this study was downloaded from UCSC Xena Browser. 
Dataset Link: https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443


## Setup
This project requires the following packages:
- 'MOFA2'
- 'mixOmics'
- clusterProfiler'
- 'survminer'
- survival
- BioCManager
