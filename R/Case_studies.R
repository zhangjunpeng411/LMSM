######## Part 1: Running scripts for identifying LncmiRSM modules, the WGCNA method is used for identifying lncRNA-mRNA coexpression modules ########
## Load required packages and utility functions
library(WGCNA)
library(PMA)
library(genefu)
library(varhandle)
library(broom)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(survival)
library(SPONGE)
library(mldr)
library(utiml)
library(e1071)
library(corpcor)
library(miRspongeR)
source("LMSM.R")

## Load data source
load("BRCA_miRNA_lncRNA_mRNA.RData")
miRTarget <- read.csv("miRTarBase_v8.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv", header = TRUE, sep = ",")
BRCA_lncRNA <- read.csv("LncRNADisease_v2.0+Lnc2Cancer_v2.0+MNDR_v2.0.csv", header = FALSE, sep = ",")
BRCA_mRNA <- read.csv("DisGeNET_v5.0+COSMIC_v86.csv", header = FALSE, sep = ",")
BRCA_gene <- rbind(BRCA_lncRNA, BRCA_mRNA)
Validated_sponge_lncRNA_experiment <- read.csv("miRSponge+LncCeRBase+LncACTdb_v2.0.csv", header = FALSE, sep = ",")
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

## Evaluate the significance of each LMSM module by using null model
LMSM_WGCNA_precomputed_cov_matrices <- precomputed_cov_matrices
LMSM_WGCNA_null_model <- sponge_build_null_model(number_of_datasets = 1e+06, number_of_samples = 500, 
                                                 cov_matrices = LMSM_WGCNA_precomputed_cov_matrices,  
						 ks = seq(0.8, 0.9, 0.1), m_max = 1)
LMSM_WGCNA_modules <- data.frame(geneA = paste("ceR_module", seq(nrow(LMSM_WGCNA_Filter_modules))), 
                                 geneB = paste("mR_module", seq(nrow(LMSM_WGCNA_Filter_modules))), 
				 df = replicate(nrow(LMSM_WGCNA_Filter_modules), 1), 
				 cor = LMSM_WGCNA_Filter_modules[, 6], 
				 pcor = LMSM_WGCNA_Filter_modules[, 9], 
				 mscor = LMSM_WGCNA_Filter_modules[, 10])
LMSM_WGCNA_modules_p.values <- sponge_compute_p_values(sponge_result = LMSM_WGCNA_modules, null_model = LMSM_WGCNA_null_model)

## BRCA subtypes identification using PAM50 method in genefu package
annot <-  data.frame(colnames(RNASeqV2_USE))
colnames(annot) <- "Gene.Symbol"
BRCA.subtype <- molecular.subtyping(sbt.model="pam50", data = RNASeqV2_USE, annot = annot, do.mapping = FALSE)$subtype
BRCA.subtype <- cbind(names(BRCA.subtype), unfactor(BRCA.subtype))
colnames(BRCA.subtype) <- c("TCGA samples", "BRCA subtype")

## Estimate GSVA enrichment score using gsva method in GSVA package
expr <- t(cbind(LncRNA_USE, RNASeqV2_USE))
BRCA_gsva_es <- gsva(expr, LMSM_WGCNA_Modulegenes, mx.diff = FALSE)
rownames(BRCA_gsva_es) <- paste("LMSM", seq_along(LMSM_WGCNA_Modulegenes), sep=" ")

## Identify BRCA subtype-specific LMSM modules
subtype_specific_module_up <- BRCA.subtype.specific.module(BRCA_gsva_es, BRCA.subtype[, 2], GSVA.enrichment.type = "Up")
subtype_specific_module_down<- BRCA.subtype.specific.module(BRCA_gsva_es, BRCA.subtype[, 2], GSVA.enrichment.type = "Down")
LMSM_subtype_specific_module_up <- cbind(rownames(LMSM_WGCNA_Filter_modules), subtype_specific_module_up)
LMSM_subtype_specific_module_down <- cbind(rownames(LMSM_WGCNA_Filter_modules), subtype_specific_module_down)
colnames(LMSM_subtype_specific_module_up) <- colnames(LMSM_subtype_specific_module_down) <- c("LMSM_modules", "subtype specific")

annotation_col <- data.frame(BRCA_subtype = BRCA.subtype[, 2])
rownames(annotation_col) <- BRCA.subtype[, 1]
BRCA_gsva_es_order <- BRCA_gsva_es[, order(BRCA.subtype[, 2])]
BRCA_gsva_es_order_up <- BRCA_gsva_es_order[which(subtype_specific_module_up != "uncertain"), ]
BRCA_gsva_es_order_down <- BRCA_gsva_es_order[which(subtype_specific_module_down != "uncertain"), ]
interin1 <- LMSM_subtype_specific_module_up[which(subtype_specific_module_up != "uncertain"), ]
interin2 <- LMSM_subtype_specific_module_down[which(subtype_specific_module_down != "uncertain"), ]
rownames(BRCA_gsva_es_order_up) <- paste(interin1[, 1], " (", interin1[, 2], "-specific)", sep="")
rownames(BRCA_gsva_es_order_down) <- paste(interin2[, 1], " (", interin2[, 2], "-specific)", sep="")
ann_colors = list(BRCA_subtype = c(Basal = "#7570B3", Her2 = "#E7298A", LumA = "#66A61E", LumB = "#1B9E77", Normal = "#D95F02"))
pheatmap(BRCA_gsva_es_order_up,  color = colorRampPalette(c("navy", "white", "firebrick3"))(300),
         annotation_col = annotation_col, annotation_colors = ann_colors, cluster_row = FALSE, cluster_col = FALSE, cutree_cols = 5, 
         show_colnames = FALSE, main = "Up-regulated BRCA subtype-specific LMSM modules")
pheatmap(BRCA_gsva_es_order_down,  color = colorRampPalette(c("navy", "white", "firebrick3"))(300),
         annotation_col = annotation_col, annotation_colors = ann_colors, cluster_row = FALSE, cluster_col = FALSE, cutree_cols = 5, 
         show_colnames = FALSE, main = "Down-regulated BRCA subtype-specific LMSM modules")

## BRCA enrichment analysis of LMSM modules
LMSM_WGCNA_Modulegenes_BRCA_EA <- module.BRCA.EA(LncRNA_USE, RNASeqV2_USE, BRCA_gene, LMSM_WGCNA_Modulegenes)

## Validation of lncRNA related miRNA sponge interactions in each LMSM module
WGCNA_validate_res_experiment <- sponge.lncRNA.validate(LMSM_WGCNA_Modulegenes, Validated_sponge_lncRNA_experiment)

## Extract lncRNA-related miRNA sponge interactions of each LMSM module
LMSM_WGCNA_miRSponge <- lapply(seq(LMSM_WGCNA_Modulegenes), function(i) Extract.miRSponge(LncRNA_USE, RNASeqV2_USE, LMSM_WGCNA_Modulegenes[[i]]))

## Extract miRNA-target interactions of each LMSM module
LMSM_WGCNA_miRTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) Extract.miRTarget(LMSM_WGCNA_CommonmiRs[[i]], LMSM_WGCNA_Modulegenes[[i]]))

## Understand predicted lncRNA-related miRNA sponge interactions, predicted and putative miRNA-target interactions of each LMSM module
LMSM_WGCNA_understand_miRSpongeTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) 
                                                    Understand.miRSpongeTarget(LncRNA_USE, RNASeqV2_USE, miRTarget, LMSM_WGCNA_CommonmiRs[[i]], 
						    LMSM_WGCNA_Modulegenes[[i]]))
LMSM_WGCNA_miRSpongeTarget <- do.call("rbind", LMSM_WGCNA_understand_miRSpongeTarget)
rownames(LMSM_WGCNA_miRSpongeTarget) <- names(LMSM_WGCNA_CommonmiRs)

## Survival analysis of each LMSM module
ExpData <- cbind(LncRNA_USE, RNASeqV2_USE)
SurvData <- data.frame(rownames(clinical_USE), clinical_USE$time, clinical_USE$status)
colnames(SurvData) <- c("sample", "time", "status")
sponge_WGCNA_Module_Survival <- moduleSurvival(LMSM_WGCNA_Modulegenes, ExpData, SurvData, devidePercentage=.5)

## miRNAs distribution in LMSM modules
LMSM_WGCNA_miR_distribution <- miR.distribution(LMSM_WGCNA_CommonmiRs)

## Performance of each LMSM module for classifying BRCA subtypes
LMSM_WGCNA_classify_baseline <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, LMSM_WGCNA_Modulegenes, method = "baseline"))
LMSM_WGCNA_classify <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, LMSM_WGCNA_Modulegenes, method = "br"))

## Comparison with a graph clustering-based strategy, the CPU runtime of the method is high. Users can divide the size of the variable "RNASeqV2_USE" into several parts (e.g. 8 parts)
SPPC_network <- SPPC(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget)
SPPC_module_genes <- netModule(SPPC_network[, 1:2], modulesize = 4)
SPPC_module <- CandModgenes(LncRNA_USE, RNASeqV2_USE, SPPC_module_genes)

# BRCA enrichment analysis
SPPC_Modulegenes_BRCA_EA <- module.BRCA.EA(LncRNA_USE, RNASeqV2_USE, BRCA_gene, SPPC_module)

# Validation analysis
SPPC_validate_res_experiment <- sponge.lncRNA.validate(SPPC_module, Validated_sponge_lncRNA_experiment)

# Survival analysis
sponge_SPPC_Module_Survival <- moduleSurvival(SPPC_module, ExpData, SurvData, devidePercentage=.5)

# Performance for classifying BRCA subtypes
SPPC_classify <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, SPPC_module, method = "br"))

save.image("LMSM_WGCNA_0.8_BRCA_miRNA_lncRNA_mRNA.RData")



######## Part 2: Running scripts for identifying LMSM modules, the SGFA method is used for identifying lncRNA-mRNA coexpression modules ########
## Load required packages and utility functions
library(GFA)
library(PMA)
library(genefu)
library(varhandle)
library(broom)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(survival)
library(SPONGE)
library(mldr)
library(utiml)
library(e1071)
library(WGCNA)
library(corpcor)
library(miRspongeR)
source("LMSM.R")

## Load data source
load("BRCA_miRNA_lncRNA_mRNA.RData")
miRTarget <- read.csv("miRTarBase_v8.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv", header = TRUE, sep = ",")
BRCA_lncRNA <- read.csv("LncRNADisease_v2.0+Lnc2Cancer_v2.0+MNDR_v2.0.csv", header = FALSE, sep = ",")
BRCA_mRNA <- read.csv("DisGeNET_v5.0+COSMIC_v86.csv", header = FALSE, sep = ",")
BRCA_gene <- rbind(BRCA_lncRNA, BRCA_mRNA)
Validated_sponge_lncRNA_experiment <- read.csv("miRSponge+LncCeRBase+LncACTdb_v2.0.csv", header = FALSE, sep = ",")
set.seed(12345)

## Idenfying lncRNA related miRNA sponge modules
CandidateModulegenes_SGFA <- module_SGFA(LncRNA_USE, RNASeqV2_USE, StrengthCut = 0.8)
CommonmiRs_SGFA <- share_miRs(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget, CandidateModulegenes_SGFA)
LMSM_SGFA <- LMSM(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget, CandidateModulegenes_SGFA)
LMSM_SGFA_Filter_modules <- LMSM_SGFA[which(LMSM_SGFA[, 3] >=3 & LMSM_SGFA[, 5] < 0.05 & 
                                              LMSM_SGFA[, 6] > 0.8 & LMSM_SGFA[, 10] > 0.1), ]
LMSM_SGFA_Modulegenes <- lapply(which(LMSM_SGFA[, 3] >=3 & LMSM_SGFA[, 5] < 0.05 & 
                                    LMSM_SGFA[, 6] > 0.8 & LMSM_SGFA[, 10] > 0.1), 
                                    function(i) CandidateModulegenes_SGFA[[i]])
LMSM_SGFA_CommonmiRs <- lapply(which(LMSM_SGFA[, 3] >=3 & LMSM_SGFA[, 5] < 0.05 & 
                                   LMSM_SGFA[, 6] > 0.8 & LMSM_SGFA[, 10] > 0.1), 
				   function(i) CommonmiRs_SGFA[[i]])
rownames(LMSM_SGFA_Filter_modules) <- names(LMSM_SGFA_Modulegenes) <- names(LMSM_SGFA_CommonmiRs) <- paste("LMSM", seq_along(LMSM_SGFA_Modulegenes), sep=" ")

## Evaluate the significance of each LMSM module by using null model
LMSM_SGFA_precomputed_cov_matrices <- precomputed_cov_matrices
LMSM_SGFA_null_model <- sponge_build_null_model(number_of_datasets = 1e+06, number_of_samples = 500, 
                                                cov_matrices = LMSM_SGFA_precomputed_cov_matrices,  
						ks = seq(0.8, 0.9, 0.1), m_max = 1)
LMSM_SGFA_modules <- data.frame(geneA = paste("ceR_module", seq(nrow(LMSM_SGFA_Filter_modules))), 
                                 geneB = paste("mR_module", seq(nrow(LMSM_SGFA_Filter_modules))), 
				 df = replicate(nrow(LMSM_SGFA_Filter_modules), 1), 
				 cor = LMSM_SGFA_Filter_modules[, 6], 
				 pcor = LMSM_SGFA_Filter_modules[, 9], 
				 mscor = LMSM_SGFA_Filter_modules[, 10])
LMSM_SGFA_modules_p.values <- sponge_compute_p_values(sponge_result = LMSM_SGFA_modules, null_model = LMSM_SGFA_null_model)

## BRCA subtypes identification using PAM50 method in genefu package
annot <-  data.frame(colnames(RNASeqV2_USE))
colnames(annot) <- "Gene.Symbol"
BRCA.subtype <- molecular.subtyping(sbt.model="pam50", data = RNASeqV2_USE, annot = annot, do.mapping = FALSE)$subtype
BRCA.subtype <- cbind(names(BRCA.subtype), unfactor(BRCA.subtype))
colnames(BRCA.subtype) <- c("TCGA samples", "BRCA subtype")

## Estimate GSVA enrichment score using gsva method in GSVA package
expr <- t(cbind(LncRNA_USE, RNASeqV2_USE))
BRCA_gsva_es <- gsva(expr, LMSM_SGFA_Modulegenes, mx.diff = FALSE)
rownames(BRCA_gsva_es) <- paste("LMSM", seq_along(LMSM_SGFA_Modulegenes), sep=" ")

## Identify BRCA subtype-specific LMSM modules
subtype_specific_module_up <- BRCA.subtype.specific.module(BRCA_gsva_es, BRCA.subtype[, 2], GSVA.enrichment.type = "Up")
subtype_specific_module_down<- BRCA.subtype.specific.module(BRCA_gsva_es, BRCA.subtype[, 2], GSVA.enrichment.type = "Down")
LMSM_subtype_specific_module_up <- cbind(rownames(LMSM_SGFA_Filter_modules), subtype_specific_module_up)
LMSM_subtype_specific_module_down <- cbind(rownames(LMSM_SGFA_Filter_modules), subtype_specific_module_down)
colnames(LMSM_subtype_specific_module_up) <- colnames(LMSM_subtype_specific_module_down) <- c("LMSM_modules", "subtype specific")

annotation_col <- data.frame(BRCA_subtype = BRCA.subtype[, 2])
rownames(annotation_col) <- BRCA.subtype[, 1]
BRCA_gsva_es_order <- BRCA_gsva_es[, order(BRCA.subtype[, 2])]
BRCA_gsva_es_order_up <- BRCA_gsva_es_order[which(subtype_specific_module_up != "uncertain"), ]
BRCA_gsva_es_order_down <- BRCA_gsva_es_order[which(subtype_specific_module_down != "uncertain"), ]
interin1 <- LMSM_subtype_specific_module_up[which(subtype_specific_module_up != "uncertain"), ]
interin2 <- LMSM_subtype_specific_module_down[which(subtype_specific_module_down != "uncertain"), ]
rownames(BRCA_gsva_es_order_up) <- paste(interin1[, 1], " (", interin1[, 2], "-specific)", sep="")
rownames(BRCA_gsva_es_order_down) <- paste(interin2[, 1], " (", interin2[, 2], "-specific)", sep="")
ann_colors = list(BRCA_subtype = c(Basal = "#7570B3", Her2 = "#E7298A", LumA = "#66A61E", LumB = "#1B9E77", Normal = "#D95F02"))
pheatmap(BRCA_gsva_es_order_up,  color = colorRampPalette(c("navy", "white", "firebrick3"))(300),
         annotation_col = annotation_col, annotation_colors = ann_colors, cluster_row = FALSE, cluster_col = FALSE, cutree_cols = 5, 
         show_colnames = FALSE, main = "Up-regulated BRCA subtype-specific LMSM modules")
pheatmap(BRCA_gsva_es_order_down,  color = colorRampPalette(c("navy", "white", "firebrick3"))(300),
         annotation_col = annotation_col, annotation_colors = ann_colors, cluster_row = FALSE, cluster_col = FALSE, cutree_cols = 5, 
         show_colnames = FALSE, main = "Down-regulated BRCA subtype-specific LMSM modules")

## BRCA enrichment analysis of LMSM modules
LMSM_SGFA_Modulegenes_BRCA_EA <- module.BRCA.EA(LncRNA_USE, RNASeqV2_USE, BRCA_gene, LMSM_SGFA_Modulegenes)

## Validation of lncRNA related miRNA sponge interactions in each LMSM module
SGFA_validate_res_experiment <- sponge.lncRNA.validate(LMSM_SGFA_Modulegenes, Validated_sponge_lncRNA_experiment)

## Extract lncRNA-related miRNA sponge interactions of each LMSM module
LMSM_SGFA_miRSponge <- lapply(seq(LMSM_SGFA_Modulegenes), function(i) Extract.miRSponge(LncRNA_USE, RNASeqV2_USE, LMSM_SGFA_Modulegenes[[i]]))

## Extract miRNA-target interactions of each LMSM module
LMSM_SGFA_miRTarget <- lapply(seq(LMSM_SGFA_CommonmiRs), function(i) Extract.miRTarget(LMSM_SGFA_CommonmiRs[[i]], LMSM_SGFA_Modulegenes[[i]]))

## Understand predicted lncRNA-related miRNA sponge interactions, predicted and putative miRNA-target interactions of each LMSM module
LMSM_SGFA_understand_miRSpongeTarget <- lapply(seq(LMSM_SGFA_CommonmiRs), function(i) 
                                                   Understand.miRSpongeTarget(LncRNA_USE, RNASeqV2_USE, miRTarget, 
						   LMSM_SGFA_CommonmiRs[[i]], LMSM_SGFA_Modulegenes[[i]]))
LMSM_SGFA_miRSpongeTarget <- do.call("rbind", LMSM_SGFA_understand_miRSpongeTarget)
rownames(LMSM_SGFA_miRSpongeTarget) <- names(LMSM_SGFA_CommonmiRs)

## Survival analysis of each LMSM module
ExpData <- cbind(LncRNA_USE, RNASeqV2_USE)
SurvData <- data.frame(rownames(clinical_USE), clinical_USE$time, clinical_USE$status)
colnames(SurvData) <- c("sample", "time", "status")
sponge_SGFA_Module_Survival <- moduleSurvival(LMSM_SGFA_Modulegenes, ExpData, SurvData, devidePercentage=.5)

## miRNAs distribution in LMSM modules
LMSM_SGFA_miR_distribution <- miR.distribution(LMSM_SGFA_CommonmiRs)

## Performance of each LMSM module for classifying BRCA subtypes
LMSM_SGFA_classify_baseline <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, LMSM_SGFA_Modulegenes, method = "baseline"))
LMSM_SGFA_classify <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, LMSM_SGFA_Modulegenes, method = "br"))

## Comparison with a graph clustering-based strategy, the CPU runtime of the method is high. Users can divide the size of the variable "RNASeqV2_USE" into several parts (e.g. 8 parts)
SPPC_network <- SPPC(miRNASeqHiseq_USE, LncRNA_USE, RNASeqV2_USE, miRTarget)
SPPC_module_genes <- netModule(SPPC_network[, 1:2], modulesize = 4)
SPPC_module <- CandModgenes(LncRNA_USE, RNASeqV2_USE, SPPC_module_genes)

# BRCA enrichment analysis
SPPC_Modulegenes_BRCA_EA <- module.BRCA.EA(LncRNA_USE, RNASeqV2_USE, BRCA_gene, SPPC_module)

# Validation analysis
SPPC_validate_res_experiment <- sponge.lncRNA.validate(SPPC_module, Validated_sponge_lncRNA_experiment)

# Survival analysis
sponge_SPPC_Module_Survival <- moduleSurvival(SPPC_module, ExpData, SurvData, devidePercentage=.5)

# Performance for classifying BRCA subtypes
SPPC_classify <- do.call(cbind, module.classify(LncRNA_USE, RNASeqV2_USE, BRCA.subtype, SPPC_module, method = "br"))

save.image("LMSM_SGFA_0.8_BRCA_miRNA_lncRNA_mRNA.RData")
