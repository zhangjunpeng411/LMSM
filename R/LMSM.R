## Identification of co-expressed gene modules using SGFA method
## ExpData: rows are samples, columns are genes.
module_SGFA <- function(ceRExp, mRExp, StrengthCut = 0.8, iter.max = 5000,
    num.ModuleceRs = 2, num.ModulemRs = 2) {

    ExpData <- list(ceRExp, mRExp)    
    names(ExpData) = c("ceRNA expression", "mRNA expression")

    # Normalize the data - here we assume that every feature is equally
    # important
    norm <- normalizeData(ExpData, type = "scaleFeatures")

    # Get the model options to detect bicluster structure
    opts <- getDefaultOpts(bicluster = TRUE)

    # Check for sampling chain convergence
    opts$convergenceCheck <- TRUE
    opts$iter.max <- iter.max

    # Infer the model
    res <- gfa(norm$train, opts = opts)

    # Extract gene index of each bicluster, using stength cutoff (absolute
    # value of association)
    BCresnum <- lapply(seq_len(dim(res$W)[2]), function(i) which(abs(res$W[,
        i]) >= StrengthCut))

    # Extract genes of each bicluster
    Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(cbind(ceRExp,
        mRExp))[BCresnum[[i]]])

    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(ceRExp))))
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(mRExp))))

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])

    return(CandidateModulegenes)
}

## Identification of co-expressed gene modules using WGCNA method
## ExpData: rows are samples, columns are genes.
module_WGCNA <- function(ceRExp, mRExp, RsquaredCut = 0.8, 
    num.ModuleceRs = 2, num.ModulemRs = 2){

    ExpData <- cbind(ceRExp, mRExp)

    Optimalpower <- pickSoftThreshold(ExpData, RsquaredCut = RsquaredCut)$powerEstimate
    adjacencymatrix <- adjacency(ExpData, power = Optimalpower)
    dissTOM <- TOMdist(adjacencymatrix)
    hierTOM <- flashClust(as.dist(dissTOM), method = "average")

    # The function cutreeDynamic colors each gene by the branches
    # that result from choosing a particular height cutoff.
    colorh <- cutreeDynamic(hierTOM, method="tree") + 1
    StandColor <- c("grey", standardColors(n = NULL))
    colorh <- unlist(lapply(seq_len(length(colorh)), function(i) StandColor[colorh[i]]))
    colorlevels <- unique(colorh)
    colorlevels <- colorlevels[-which(colorlevels=="grey")]

    Modulegenes <- lapply(seq_len(length(colorlevels)), function(i)
                          colnames(ExpData)[ which(colorh==colorlevels[i]) ])

    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in% colnames(ceRExp))) )
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in% colnames(mRExp))) )

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])

    return(CandidateModulegenes)
}

## Extract common miRNAs for each co-expressed gene module
share_miRs <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes){

        miRTarget <- as.matrix(miRTarget)            
        miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% c(colnames(ceRExp),colnames(mRExp)))),]

	Res <- list()

        for (i in seq_along(CandidateModulegenes)){
            
            tmp1 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
	                   intersect(CandidateModulegenes[[i]], colnames(ceRExp)) ), 1])        
            
            tmp2 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
	                   intersect(CandidateModulegenes[[i]], colnames(mRExp)) ), 1])
            
            tmp3 <- intersect( tmp1, tmp2 )

	    Res[[i]] <- tmp3
	    
	}

return(Res)
}

## Evaluate candidate ceRM-mRM pairs as miRNA sponge modules using sensitivity canonical correlation
LMSM <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes){   
         
        miRTarget <- as.matrix(miRTarget)            
        miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% c(colnames(ceRExp),colnames(mRExp)))),]
        Res <- c()

        for (i in seq_along(CandidateModulegenes)){
            # Calculate significance of miRNAs shared by each ceRM-mRM pair
            tmp1 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
	                   intersect(CandidateModulegenes[[i]], colnames(ceRExp)) ), 1])        
            M1 <- length(tmp1)
            tmp2 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
	                   intersect(CandidateModulegenes[[i]], colnames(mRExp)) ), 1])
            M2 <- length(tmp2)
            tmp3 <- intersect( tmp1, tmp2 )
            M3 <- length( tmp3 )
            M4 <- length( colnames(miRExp) )        
            M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

            if(M3>=3){                    
            
            # Canonical correlation between a group of ceRNAs and a group of mRNAs
            perm.out_ceR_mR <- CCA.permute(ceRExp[, which(colnames(ceRExp) %in% 
	                                   CandidateModulegenes[[i]])], mRExp[, which(colnames(mRExp) %in% 
					   CandidateModulegenes[[i]])], typex="standard", typez="standard", nperms=100)
            out_ceR_mR <- CCA(ceRExp[, which(colnames(ceRExp) %in% 
	                      CandidateModulegenes[[i]])], mRExp[, which(colnames(mRExp) %in% 
			      CandidateModulegenes[[i]])], typex="standard", typez="standard", 
			      K=1, penaltyx=perm.out_ceR_mR$bestpenaltyx, 
			      penaltyz=perm.out_ceR_mR$bestpenaltyz, v=perm.out_ceR_mR$v.init)
	    M6 <- out_ceR_mR$cor
 
            # Canonical correlation between a group of miRNAs and a group of mRNAs
            perm.out_miR_mR <- CCA.permute(miRExp[, which(colnames(miRExp) %in% tmp3)], 
	                                   mRExp[, which(colnames(mRExp) %in% 
					   CandidateModulegenes[[i]])], typex="standard", typez="standard", nperms=100)
            out_miR_mR <- CCA(miRExp[, which(colnames(miRExp) %in% tmp3)], 
	                  mRExp[, which(colnames(mRExp) %in% CandidateModulegenes[[i]])], 
			  typex="standard", typez="standard", K=1, penaltyx=perm.out_miR_mR$bestpenaltyx, 
			  penaltyz=perm.out_miR_mR$bestpenaltyz, v=perm.out_miR_mR$v.init)	    
            M7 <- out_miR_mR$cor      
 
            # Canonical correlation between a group of miRNAs and a group of ceRNAs
            perm.out_miR_ceR <- CCA.permute(miRExp[, which(colnames(miRExp) %in% tmp3)], 
	                                    ceRExp[, which(colnames(ceRExp) %in% CandidateModulegenes[[i]])], 
	                                    typex="standard", typez="standard", nperms=100)
            out_miR_ceR <- CCA(miRExp[, which(colnames(miRExp) %in% tmp3)], 
	                       ceRExp[, which(colnames(ceRExp) %in% CandidateModulegenes[[i]])], 
			       typex="standard", typez="standard", K=1, 
			       penaltyx=perm.out_miR_ceR$bestpenaltyx, 
			       penaltyz=perm.out_miR_ceR$bestpenaltyz, v=perm.out_miR_ceR$v.init)
	    M8 <- out_miR_ceR$cor

            # Calculate partial canonical correlation between a group of ceRNAs and a group of mRNAs on condition a group of miRNAs
            M9 <- (M6 - M7*M8)/(sqrt(1 - M7^2)*sqrt(1 - M8^2))

            # Calculate sensitivity canonical correlation between a group of ceRNAs and a group of mRNAs on condition a group of miRNAs
            M10 <- M6 - M9 
           
            } else {

            M6 <- NA; M7 <- NA; M8 <- NA; M9 <- NA; M10 <- NA   
        
            }

            tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)    
            Res <- rbind(Res, tmp)
        }
            colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs", "#Shared miRNAs", 
                               "#Background miRNAs", "Sig. p.value of sharing miRNAs", "Canonical correlation of ceRM-mRM pairs", 
                               "Canonical correlation of miRM-mRM pairs", "Canonical correlation of miRM-ceRM pairs", 
                               "Partial canonical correlation of ceRM-mRM pairs", "Sensitivity canonical correlation of ceRM-mRM pairs")
            rownames(Res) <- paste("Module", seq_along(CandidateModulegenes), sep=" ")
return(Res)
}

## Determine the LMSM modules belong to specific BRCA subtype
BRCA.subtype.specific.module <- function(BRCA_gsva_es, BRCA.subtype, 
    p.adjust.method = "BH", p.value.cutoff = 0.05, GSVA.enrichment.type = c("Up", "Down")){

    subtype_specific_module <- c()
    for (i in seq_len(dim(BRCA_gsva_es)[1])){
        Response <- BRCA_gsva_es[i, ]
        Treatment <- BRCA.subtype    
        RT <- tidy(pairwise.t.test(Response, Treatment, p.adjust.method = p.adjust.method))
        BRCA.subtype.unique <- unique(Treatment)
        mu <- unlist(lapply(seq(BRCA.subtype.unique), function(i) mean(Response[Treatment==BRCA.subtype.unique[i]]))) 
        if (GSVA.enrichment.type == "Up") {  
        BRCA_subtype_specific <- BRCA.subtype.unique[which(mu == max(mu))]
        } else if (GSVA.enrichment.type == "Down") {
        BRCA_subtype_specific <- BRCA.subtype.unique[which(mu == min(mu))]
        }
        BRCA_subtype_specific_p.value <- RT$p.value[c(which(RT$group1 == BRCA_subtype_specific), 
                                             which(RT$group2 == BRCA_subtype_specific))]
        if(all(BRCA_subtype_specific_p.value < p.value.cutoff)){
            subtype_specific_module[i] <- BRCA_subtype_specific
        } else {
            subtype_specific_module[i] <- "uncertain"
        }
    }
return(subtype_specific_module)
} 

## BRCA enrichment analysis using hypergeometric distribution test 
module.BRCA.EA <- function(ceRExp, mRExp, BRCAgenes, Modulelist) {

    ExpData <- cbind(ceRExp, mRExp)      

    B <- ncol(ExpData)
    N <- length(intersect(colnames(ExpData), as.matrix(BRCAgenes)))
    M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
    x <- unlist(lapply(seq_along(Modulelist), function(i) length(intersect(Modulelist[[i]], as.matrix(BRCAgenes)))))    
    p.value <- 1 - phyper(x - 1, N, B - N, M)
    
    names(p.value) <- names(Modulelist)
    return(p.value)
}

## Validation of lncRNA related miRNA sponge interactions in each LMSM module
sponge.lncRNA.validate <- function(Modulelist, Validated_sponge_lncRNA) {

    validate_res <- lapply(seq(Modulelist), function(i) Validated_sponge_lncRNA[intersect(which(as.matrix(Validated_sponge_lncRNA[, 1]) %in% Modulelist[[i]]), 
                        which(as.matrix(Validated_sponge_lncRNA[, 2]) %in% Modulelist[[i]])), ])
    return(validate_res)
}

## Extract lncRNA-related miRNA sponge interactions of each LMSM module
Extract.miRSponge <- function(ceRExp, mRExp, Modulegenes){
    
    ceRNAs <- Modulegenes[which(Modulegenes %in% colnames(ceRExp))]
    mRNAs <- Modulegenes[which(Modulegenes %in% colnames(mRExp))]
    len_ceRNAs <- length(ceRNAs)
    len_mRNAs <- length(mRNAs)
    res_int <- matrix(NA, len_ceRNAs*len_mRNAs, 2)
    for (i in seq_len(len_ceRNAs)){
        for (j in seq_len(len_mRNAs)){
	    res_int[(i-1)*len_mRNAs+j, 1] <- ceRNAs[i]
            res_int[(i-1)*len_mRNAs+j, 2] <- mRNAs[j]
	}
    }
return(res_int)
}

## Extract miRNA-target interactions of each LMSM module
Extract.miRTarget <- function(CommonmiRs, targets){

    len_CommonmiRs <- length(CommonmiRs)
    len_targets <- length(targets)
    res_int <- matrix(NA, len_CommonmiRs*len_targets, 2)
    for (i in seq_len(len_CommonmiRs)){
        for (j in seq_len(len_targets)){
	    res_int[(i-1)*len_targets+j, 1] <- CommonmiRs[i]
            res_int[(i-1)*len_targets+j, 2] <- targets[j]
	}
    }
return(res_int)
}

## Understand predicted lncRNA-related miRNA sponge interactions, predicted and putative miRNA-target interactions of each LMSM module
Understand.miRSpongeTarget <- function(ceRExp, mRExp, miRTarget, CommonmiRs, targets){

    len_CommonmiRs <- length(CommonmiRs)
    lncRNAs <- colnames(ceRExp)[which(colnames(ceRExp) %in% targets)]
    len_lncRNAs <- length(lncRNAs)
    mRNAs <- colnames(mRExp)[which(colnames(mRExp) %in% targets)]
    len_mRNAs <- length(mRNAs)
    num_ceRmR <- len_lncRNAs * len_mRNAs
    num_miRlncR <- len_CommonmiRs * len_lncRNAs
    num_miRmR <- len_CommonmiRs * len_mRNAs
    num_putative_miRlncR <- dim(miRTarget[intersect(which(miRTarget[, 1] %in% CommonmiRs),
	                      which(miRTarget[, 2] %in% lncRNAs)),])[1]
    num_putative_miRmR <- dim(miRTarget[intersect(which(miRTarget[, 1] %in% CommonmiRs),
	                      which(miRTarget[, 2] %in% mRNAs)),])[1]
    res <- c(len_CommonmiRs, len_lncRNAs, len_mRNAs, num_ceRmR, num_miRlncR, num_miRmR, num_putative_miRlncR,  num_putative_miRmR)
    names(res) <- c("#Shared miRNAs", "#lncRNAs", "#mRNAs", "#Predicted lncRNA-related miRNA sponge interactions", "#Predicted miRNA-lncRNA", "#Predicted miRNA-mRNA", "#Putative miRNA-lncRNA", "#Putative miRNA-mRNA")

return(res)
}

# Survival analysis of each LMSM module
moduleSurvival <- function (Modulelist, ExpData, SurvData, devidePercentage = 0.5, 
    plot = FALSE) 
{
    ExpDataNames <- colnames(ExpData)
    myfit <- list()
    LogRank <- list()
    for (i in seq_along(Modulelist)) {
        Interin_Data <- cbind(SurvData[, seq(2, 3)], ExpData[, 
            which(ExpDataNames %in% Modulelist[[i]])])
        Interin_Data <- na.omit(Interin_Data)
        try_mm <- try(coxph(survival::Surv(time, status) ~ ., 
            data = data.frame(Interin_Data)), silent = TRUE)
        if ("try-error" %in% class(try_mm)) 
            next
        mm <- coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data))
        Risk_score <- predict(mm, newdata = data.frame(Interin_Data), 
            type = "risk")
        group <- rep("NA", dim(Interin_Data)[1])
        group[Risk_score > quantile(Risk_score, probs = devidePercentage)] <- "High"
        group[Risk_score <= quantile(Risk_score, probs = devidePercentage)] <- "Low"
        Data <- cbind(Interin_Data[, seq_len(2)], group)
        myfit[[i]] <- survfit(survival::Surv(time, status) ~ 
            group, data = Data)
        sdf <- survdiff(survival::Surv(time, status) ~ group, 
            data = Data)
        sdf.p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        HR <- (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
        HRlow95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/sdf$exp[1] + 
            1/sdf$exp[2]))
        HRup95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/sdf$exp[1] + 
            1/sdf$exp[2]))
        LogRank[[i]] <- c(sdf$chisq, sdf.p.val, HR, HRlow95, 
            HRup95)
    }
    if (plot) {
        for (i in seq_along(myfit)) {
            if (!is.null(LogRank[[i]])) {
                dev.new()
                plot(myfit[[i]], lty = 1, col = c("red", 
                  "green"), main = paste("Module", 
                  i), xlab = "Time (Months)", ylab = "Probability of survival")
                legend("topright", legend = c("High risk group", 
                  "Low risk group"), lty = seq_len(2), 
                  col = c("red", "green"))
            }
        }
    }
    LogRank_res <- do.call(rbind, LogRank)
    if (length(myfit) >= 1) {
        colnames(LogRank_res) <- c("Chi-square", "p-value", 
            "HR", "HRlow95", "HRup95")
        names(LogRank) <- seq_along(Modulelist)
        LogRank[sapply(LogRank, is.null)] <- NULL
        rownames(LogRank_res) <- paste("Module", names(LogRank))
    }
    return(LogRank_res)
}

## miRNAs distribution in LMSM modules
miR.distribution <- function(CommonmiRslist) {

    miRs <- unique(unlist(CommonmiRslist))
    res <- NULL
    interin <- NULL
    for (i in seq_along(miRs)) {
        for (j in seq_along(CommonmiRslist)) {
	    if (length(which(miRs[i] %in% CommonmiRslist[[j]]) == 1)) {
	        interin <- c(interin, names(CommonmiRslist)[j])
	    }
	}
	res1 <- paste(interin, collapse = ", ")        
	res2 <- length(interin)
        res <- rbind(res, c(miRs[i], res1, res2))
        interin <- NULL
    }

return(res)
}

## Evaluate the performance of each LMSM module for classifying BRCA subtypes
module.classify <- function(ceRExp, mRExp, BRCA_subtype, Modulelist, method = "br", base.algorith = "SVM", cv.folds = 10, 
	                    cv.sampling = "stratified", cv.seed = 12345) {

    module_ceRExp <- lapply(seq_along(Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]])])
    module_mRExp <- lapply(seq_along(Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist[[i]])])
    Basal <- as.numeric(BRCA_subtype[, 2] == "Basal")
    Her2 <- as.numeric(BRCA_subtype[, 2] == "Her2")
    LumA <- as.numeric(BRCA_subtype[, 2] == "LumA")
    LumB <- as.numeric(BRCA_subtype[, 2] == "LumB")
    Normal <- as.numeric(BRCA_subtype[, 2] == "Normal")    
    module_classify <- list()

    for (i in seq_along(Modulelist)){        
	
        temp <- as.data.frame(cbind(module_ceRExp[[i]], module_mRExp[[i]], Basal, Her2, LumA, LumB, Normal))
	Indices <- ncol(temp)
	temp_mldr <- mldr_from_dataframe(temp, labelIndices = c(Indices-4, Indices-3, Indices-2, Indices-1, Indices), name = "TEMPMLDR")
        temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds, 
	               cv.sampling = cv.sampling, cv.seed = cv.seed)
        module_classify[[i]] <- temp_res

    }

    return(module_classify)
}

# Identifying lncRNA related miRNA sponge interactions using Sensitivity Partial Pearson Correlation (SPPC) method
SPPC <- function(miRExp, ceRExp, mRExp, miRTarget, minSharedmiR = 3, pvaluecutoff = 0.05, poscorcutoff = 0, senscorcutoff = 0.1){   
         
        miRTarget <- as.matrix(miRTarget)	
        miRceR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(ceRExp))), ]
        miRmR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(mRExp))), ]

	ceRSym <- unique(miRceR[, 2])
        mRSym <- unique(miRmR[, 2])
	miRSym <- unique(c(miRceR[, 1], miRmR[, 1]))

	m <- length(ceRSym)
	n <- length(mRSym)

        ceRExp_query <- ceRExp[, which(colnames(ceRExp) %in% ceRSym)]
        mRExp_query <- mRExp[, which(colnames(mRExp) %in% mRSym)]
	
        Cor.Pvalue <- corAndPvalue(ceRExp_query, mRExp_query)
        
	index <- which(Cor.Pvalue$cor > poscorcutoff & Cor.Pvalue$p < pvaluecutoff, arr.ind = TRUE)
        
	# Initialize variables
        Res <- c()

        for (i in seq_len(nrow(index))) {
	    
	        Interin1 <- miRceR[which(miRceR[, 2] %in% ceRSym[index[i, 1]]), 1]
                Interin2 <- miRmR[which(miRmR[, 2] %in% mRSym[index[i, 2]]), 1]

		M1 <- length(Interin1)
                M2 <- length(Interin2)
                M3 <- length(intersect(Interin1, Interin2))
                M4 <- length(miRSym)
                M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

                if (M3 >= minSharedmiR & M5 < pvaluecutoff) {
                    
                    C1 <- ceRSym[index[i, 1]]
                    C2 <- mRSym[index[i, 2]]

                    ceRExpIdx <- which(colnames(ceRExp) %in% ceRSym[index[i, 1]])
                    mRExpIdx <- which(colnames(mRExp) %in% mRSym[index[i, 2]])
                    miRExpIdx <- which(colnames(miRExp) %in% intersect(Interin1, Interin2))

		    M6 <- Cor.Pvalue$cor[index[i, 1], index[i, 2]]
		    M7 <- Cor.Pvalue$p[index[i, 1], index[i, 2]]
                    M8 <- M6 - corpcor::pcor.shrink(cbind(ceRExp[, ceRExpIdx], mRExp[, mRExpIdx],
                                                    miRExp[, miRExpIdx]), verbose = FALSE)[1, 2]
                } else {

                C1 <- NA; C2 <- NA; M6 <- NA; M7 <- NA; M8 <- NA 
	       
                }
	       
	        tmp <- c(C1, C2, M3, M5, M6, M7, M8)    
                Res <- rbind(Res, tmp)
	}

        Res <- Res[which(as.numeric(Res[, 7]) > senscorcutoff), ]

        colnames(Res) <- c("sponge_1", "sponge_2", "#Shared miRNAs", 
                           "Sig. p.value of sharing miRNAs", 
			   "Correlation", 
                           "Sig. p.value of correlation", 			   
			   "Sensitivity correlation")
			   
       rownames(Res) <- seq_len(nrow(Res))

return(Res)

}

## Extracting modules with at least num.ModuleceRs (e.g. 2) ceRNAs and num.ModulemRs (e.g. 2) mRNAs
CandModgenes <- function(ceRExp, mRExp, Modulegenes, num.ModuleceRs = 2, num.ModulemRs = 2){
  
    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(ceRExp))))
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(mRExp))))

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])
    
return(CandidateModulegenes)
}

