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
source("LMSM.R")

## Load data source
load("BRCA_miRNA_lncRNA_mRNA.RData")
miRTarget <- read.csv("miRTarBase_v7.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv", header = TRUE, sep = ",")
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

## Calculate average correlation between sponge lncRNAs and mRNAs in LMSM modules
LMSM_WGCNA_mean_cor <- module.avg.cor(LncRNA_USE, RNASeqV2_USE, LMSM_WGCNA_Modulegenes, resample = 1000, method = "mean")
LMSM_WGCNA_median_cor <- module.avg.cor(LncRNA_USE, RNASeqV2_USE, LMSM_WGCNA_Modulegenes, resample = 1000, method = "median")

LMSM_Random_module_mean_compare <- t.test(LMSM_WGCNA_mean_cor[[1]], LMSM_WGCNA_mean_cor[[2]], paired = TRUE)
LMSM_Random_module_median_compare <- t.test(LMSM_WGCNA_median_cor[[1]], LMSM_WGCNA_median_cor[[2]], paired = TRUE)

data1 <- data.frame(cor=c(LMSM_WGCNA_mean_cor[[1]], LMSM_WGCNA_mean_cor[[2]]), type = rep(c("LMSM modules", "Random modules"), each = 17))

p1 <- ggplot(data1)+geom_point(aes(x=rep(seq_len(17),2),y=cor,colour=type), size=3)+
          labs(x="Module ID", y="Mean value of absolute correlation", title = "p-value=1.53E-05")+
	  theme(plot.title = element_text(hjust = 0.5), legend.position=c(0.2,0.9), 
	  axis.title = element_text(face="bold"), title = element_text(face="bold"), 
          axis.text = element_text(face="bold", color = "black"), legend.text =element_text(face="bold"))

data2 <- data.frame(cor=c(LMSM_WGCNA_median_cor[[1]], LMSM_WGCNA_median_cor[[2]]), type = rep(c("LMSM modules", "Random modules"), each = 17))

p2 <- ggplot(data2)+geom_point(aes(x=rep(seq_len(17),2),y=cor,colour=type), size=3)+
          labs(x="Module ID", y="Median value of absolute correlation", title = "p-value=1.19E-05")+
	  theme(plot.title = element_text(hjust = 0.5), legend.position=c(0.2,0.9), 
	  axis.title = element_text(face="bold"), title = element_text(face="bold"), 
          axis.text = element_text(face="bold", color = "black"), legend.text =element_text(face="bold"))

## miRNAs distribution in LMSM modules
LMSM_WGCNA_miR_distribution <- miR.distribution(LMSM_WGCNA_CommonmiRs)

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
source("LMSM.R")

## Load data source
load("BRCA_miRNA_lncRNA_mRNA.RData")
miRTarget <- read.csv("miRTarBase_v7.0+TarBase_v7.0+miRWalk_v2.0+NPInter_v3.0+LncBase_v2.csv", header = TRUE, sep = ",")
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

## Calculate average correlation between sponge lncRNAs and mRNAs in LMSM modules
LMSM_SGFA_mean_cor <- module.avg.cor(LncRNA_USE, RNASeqV2_USE, LMSM_SGFA_Modulegenes, resample = 1000, method = "mean")
LMSM_SGFA_median_cor <- module.avg.cor(LncRNA_USE, RNASeqV2_USE, LMSM_SGFA_Modulegenes, resample = 1000, method = "median")

LMSM_Random_module_mean_compare <- t.test(LMSM_SGFA_mean_cor[[1]], LMSM_SGFA_mean_cor[[2]], paired = TRUE)
LMSM_Random_module_median_compare <- t.test(LMSM_SGFA_median_cor[[1]], LMSM_SGFA_median_cor[[2]], paired = TRUE)

data1 <- data.frame(cor=c(LMSM_SGFA_mean_cor[[1]], LMSM_SGFA_mean_cor[[2]]), type = rep(c("LMSM modules", "Random modules"), each = 51))

p1 <- ggplot(data1)+geom_point(aes(x=rep(seq_len(51),2),y=cor,colour=type), size=3)+
          labs(x="Module ID", y="Mean value of absolute correlation", title = "p-value=1.53E-14")+
	  theme(plot.title = element_text(hjust = 0.5), legend.position=c(0.2,0.9), 
	  axis.title = element_text(face="bold"), title = element_text(face="bold"), 
          axis.text = element_text(face="bold", color = "black"), legend.text =element_text(face="bold"))+
	  scale_y_continuous(breaks=seq(0.1, 0.3, 0.05))

data2 <- data.frame(cor=c(LMSM_SGFA_median_cor[[1]], LMSM_SGFA_median_cor[[2]]), type = rep(c("LMSM modules", "Random modules"), each = 51))

p2 <- ggplot(data2)+geom_point(aes(x=rep(seq_len(51),2),y=cor,colour=type), size=3)+
             labs(x="Module ID", y="Median value of absolute correlation", title = "p-value=5.28E-14")+
	     theme(plot.title = element_text(hjust = 0.5), legend.position=c(0.2,0.9), 
	     axis.title = element_text(face="bold"), title = element_text(face="bold"), 
             axis.text = element_text(face="bold", color = "black"), legend.text =element_text(face="bold"))+
	     scale_y_continuous(breaks=seq(0.1, 0.3, 0.05))

## miRNAs distribution in LMSM modules
LMSM_SGFA_miR_distribution <- miR.distribution(LMSM_SGFA_CommonmiRs)

save.image("LMSM_SGFA_0.8_BRCA_miRNA_lncRNA_mRNA.RData")
