# LMSM
A modular approach to identify lncRNA related miRNA sponge modules in breast cancer

## Background
Until now, existing methods for identifying lncRNA related miRNA sponge modules mainly rely on lncRNA related miRNA sponge interaction networks, which may not provide a full picture of miRNA sponging activities in biological conditions. Hence there is a strong need of new computational methods to identify lncRNA related miRNA sponge modules. In this work, we propose a framework, LMSM, to identify LncRNA related MiRNA Sponge Modules from heterogeneous data. To understand the miRNA sponging activities in biological conditions, LMSM uses gene expression data to evaluate the influence of the shared miRNAs on the clustered sponge lncRNAs and mRNAs.

## Description of each file
BRCA_miRNA_lncRNA_mRNA.RData: Matched miRNA, lncRNA and mRNA expression data, and clinical information in BRCA.

miRTarBase_v8.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv: Putative miRNA-target interactions.

LncRNADisease_v2.0+Lnc2Cancer_v2.0+MNDR_v2.0.csv: BRCA-related lncRNAs.

DisGeNET_v5.0+COSMIC_v86.csv: BRCA-related mRNAs.

miRSponge+LncCeRBase+LncACTdb_v2.0.csv: Experimentally validated lncRNA-related miRNA sponge interactions. 

LMSM.R: Functions for identifying and analyzing lncRNA related miRNA sponge modules.

Case_studies.R: Scripts of two case studies for identifying and analyzing lncRNA related miRNA sponge modules.

## Case studies
Paste all files including scripts and datasets into a single folder (set the folder as the directory of R environment), the scripts of two case studies using LMSM is implemented in Case_studies.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("Case_studies.R")
```

## The usage of LMSM
For identifying lncRNA related miRNA sponge modules, users should prepare datasets including matched miRNA, lncRNA and mRNA expression data, and putative miRNA-target (miRNA-lncRNA and miRNA-mRNA) interactions. Paste the datasets and our source file (LMSM.R) into a single folder (set the folder as the directory of R environment), users can use the following scripts to identify lncRNA related miRNA sponge modules. For convenience, the datasets prepared for users are from our datasets (BRCA_miRNA_lncRNA_mRNA.RData and miRTarBase_v8.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv). 

```{r echo=FALSE, results='hide', message=FALSE}
## Load required packages and utility functions
library(WGCNA)
library(PMA)
source("LMSM.R")

## Load data source
load("BRCA_miRNA_lncRNA_mRNA.RData")
miRTarget <- read.csv("miRTarBase_v8.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv", header = TRUE, sep = ",")
set.seed(12345)

## Idenfying lncRNA related miRNA sponge modules
CandidateModulegenes_WGCNA <- module_WGCNA(LncRNA_USE, RNASeqV2_USE, RsquaredCut = 0.8)
CommonmiRs_WGCNA <- share_miRs(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget, CandidateModulegenes_WGCNA)
LMSM_WGCNA <- LMSM(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget, CandidateModulegenes_WGCNA)
LMSM_WGCNA_Filter_modules <- LMSM_WGCNA[which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                       LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.1), ]
LMSM_WGCNA_Modulegenes <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.1), 
                                function(i) CandidateModulegenes_WGCNA[[i]])
LMSM_WGCNA_CommonmiRs <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                               LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.1), 
				                       function(i) CommonmiRs_WGCNA[[i]])
rownames(LMSM_WGCNA_Filter_modules) <- names(LMSM_WGCNA_Modulegenes) <- names(LMSM_WGCNA_CommonmiRs) <- paste("LMSM", seq_along(LMSM_WGCNA_Modulegenes), sep=" ")
```
