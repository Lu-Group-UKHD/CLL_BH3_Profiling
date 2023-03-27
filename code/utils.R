#Utility functions
library(latex2exp)

#set the global ggplot theme
theme_full <- theme_bw() + theme(axis.text = element_text(size=18),
                             axis.title = element_text(size=18),
                             axis.line = element_blank(),
                             panel.border = element_rect(size=1.5),
                             axis.ticks = element_line(size=1.5),
                             plot.title = element_text(size = 20, hjust =0.5, face="bold"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

theme_half <- theme_bw() + theme(axis.text = element_text(size=15),
                                 axis.title = element_text(size=16),
                                 axis.line =element_line(size=0.8),
                                 panel.border = element_blank(),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())

# Defien a color scheme, based on ggsci_NEJM panel, for the paper
colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")

# Function for scale a data matrix
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    x.scaled <- t(x)
  }

  if (!is.null(censor)) {
    if (length(censor) == 1) {
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      x.scaled[x.scaled > censor[2]] <- censor[2] #higher limit
      x.scaled[x.scaled < censor[1]] <- censor[1] #lower limit
    }
  }
  return(t(as.matrix(x.scaled)))
}
plotCorScatter <- function(inputTab, x, y, x_lab = "X", y_lab = "Y", title = "",
                           col = NULL, showR2 = TRUE, annoPos = "right",
                           dotCol = colList, textCol="darkred") {

  #prepare table for plotting
  plotTab <- tibble(x = inputTab[[x]],y=inputTab[[y]])
  if (!is.null(col)) plotTab <- mutate(plotTab, status = inputTab[[col]])
  plotTab <- filter(plotTab, !is.na(x), !is.na(y))

  #prepare annotation values
  corRes <- cor.test(plotTab$x, plotTab$y)
  pval <- formatNum(corRes$p.value, digits = 1, format = "e")
  Rval <- formatNum(corRes$estimate, digits = 1, format = "e")
  R2val <- formatNum(corRes$estimate^2, digits = 1, format = "e")
  Nval <- nrow(plotTab)
  annoP <- bquote(italic("P")~"="~.(pval))

  if (showR2) {
    annoCoef <-  bquote(R^2~"="~.(R2val))
  } else {
    annoCoef <- bquote(R~"="~.(Rval))
  }
  annoN <- bquote(N~"="~.(Nval))

  corPlot <- ggplot(plotTab, aes(x = x, y = y))

  if (!is.null(col)) {
    corPlot <- corPlot + geom_point(aes(fill = status), shape =21, size =3) +
      scale_fill_manual(values = dotCol)
  } else {
    corPlot <- corPlot + geom_point(fill = dotCol[1], shape =21, size=3)
  }

  corPlot <- corPlot +   geom_smooth(formula = y~x,method = "lm", se=FALSE, color = "grey50", linetype ="dashed" )

  if (annoPos == "right") {

    corPlot <- corPlot + annotate("text", x = max(plotTab$x), y = Inf, label = annoN,
                                  hjust=1, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoP,
               hjust=1, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoCoef,
               hjust=1, vjust =6, size = 5, parse = FALSE, col= textCol)

  } else if (annoPos== "left") {
    corPlot <- corPlot + annotate("text", x = min(plotTab$x), y = Inf, label = annoN,
                                  hjust=0, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoP,
               hjust=0, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoCoef,
               hjust=0, vjust =6, size = 5, parse = FALSE, col= textCol)
  }
  corPlot <- corPlot + ylab(y_lab) + xlab(x_lab) + ggtitle(title) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_full
  corPlot
}


#function for cox regression
com <- function(response, time, endpoint, scale =FALSE) {

  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)


  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}

# Function for Kaplan-Meier plot
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (years)",
               table_ratio = c(0.7,0.3), yLabelAdjust = 0) {
  #function for km plot
  survS <- tibble(time = time,
                  endpoint = endpoint)

  if (!is.null(maxTime))
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))

  if (stat == "maxstat") {
    ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                       data = survS,
                       smethod = "LogRank",
                       minprop = 0.2,
                       maxprop = 0.8,
                       alpha = NULL)

    survS$group <- factor(ifelse(response >= ms$estimate, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "median") {
    med <- median(response, na.rm = TRUE)
    survS$group <- factor(ifelse(response >= med, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "binary") {
    survS$group <- factor(response)
    if (nlevels(survS$group) > 2) {
      sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
      p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    } else {
      p <- com(survS$group, survS$time, survS$endpoint)$p
    }
  }

  if (is.null(pval)) {
    if(p< 1e-16) {
      pAnno <- bquote(italic("P")~"< 1e-16")
    } else {
      pval <- formatNum(p, digits = 1)
      pAnno <- bquote(italic("P")~"="~.(pval))
    }

  } else {
     pval <- formatNum(pval, digits = 1)
     pAnno <- bquote(italic("P")~"="~.(pval))
  }

  if (!showP) pAnno <- ""

  colorPal <- colList[1:length(unique(survS$group))]
  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE, palette = colorPal,
                  legend = ifelse(showTable, "none","top"),
                  ylab = "Fraction", xlab = "Time (years)", title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, legend.labs = sort(unique(survS$group)),
                  ggtheme = theme_half + theme(plot.title = element_text(hjust =0.5),
                                               panel.border = element_blank(),
                                               axis.title.y = element_text(vjust =yLabelAdjust)))
  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5)
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size=5)
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}

#function for plot hazard ratio
plotHazard <- function(survRes, title = "") {
  sumTab <- summary(survRes)$coefficients
  confTab <- summary(survRes)$conf.int
  #correct feature name
  nameOri <- rownames(sumTab)
  nameMod <- substr(nameOri, 1, nchar(nameOri) -1)
  plotTab <- tibble(feature = rownames(sumTab),
                    nameMod = substr(nameOri, 1, nchar(nameOri) -1),
                    HR = sumTab[,2],
                    p = sumTab[,5],
                    Upper = confTab[,4],
                    Lower = confTab[,3]) %>%
    mutate(feature = ifelse(nameMod %in% names(survRes$xlevels), nameMod, feature)) %>%
    mutate(feature = str_replace(feature, "[.]","/")) %>%
    mutate(feature = str_replace(feature, "[_]","-")) %>%
    arrange(desc(abs(p))) %>% mutate(feature = factor(feature, levels = feature)) %>%
    mutate(type = ifelse(HR >1 ,"up","down")) %>%
    mutate(Upper = ifelse(Upper > 10, 10, Upper))

  ggplot(plotTab, aes(x=feature, y = HR, color = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    geom_point(position = position_dodge(width=0.8), size=3, color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1,color = "grey20") +
    geom_text(position = position_nudge(x = 0.3),
              aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                           formatNum(p, digits = 1))),
              color = "black", size =5, parse = TRUE) +
    expand_limits(y=c(-0.5,0))+
    scale_color_manual(values = c(up = colList[1], down = colList[2])) +
    ggtitle(title) + scale_y_log10() +
    ylab("Hazard ratio") +
    coord_flip() +
    theme_full +
    theme(legend.position = "none", axis.title.y = element_blank())
}

#Function to run multivariate Cox model on test table
runCox <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Function to calculate C-index
calcConcord <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- concordance(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Function to generate pretty scientific notation format for plot label
sciPretty <- function(n, digits = 2, bold = FALSE) {
  nForm <- strsplit(format(n, digits = digits, scientific = TRUE),split = "e")
  b <- nForm[[1]][1]
  i <- as.integer(nForm[[1]][2])
  #bquote(.(b)%*%10^.(i))
  if(bold) {
    sprintf("bold(%s%%*%%10^%s)",b,i)
  } else sprintf("%s%%*%%10^%s",b,i)
}

#function to remove highly correlated values and keep track of the removing
removeCorrelated <- function(x, cutoff = 0.6, method = "pearson", keep = NULL, record = TRUE) {

  if (!is.null(keep)) {
    #if specified some feature to keep, then reorder the matrix, to make sure they are not collapsed
    posNow <- grepl(paste(keep, collapse = "|"), colnames(x))
    posNew <- rev(c(colnames(x)[posNow],colnames(x)[!posNow]))
    x <- x[,posNew]
  }
  #input is a feature matrix, features are in columns
  if (method == "binary") {
    #use binary similarity if input is a binary matrix,
    #maybe also usefull is the input is a sparse matrix
    simiMat <- 1 - as.matrix(dist(t(x), method = "binary"))
  } else if (method == "pearson") {
    #otherwise, using pearson correlation
    simiMat <- cor(x)
  } else if (method == "euclidean") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "euclidean"))
  } else if (method == "cosine") {
    # cosine similarity maybe prefered for sparse matrix
    cosineSimi <- function(x){
      x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    simiMat <- cosineSimi(t(x))
  } else if (method == "canberra") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "canberra"))/nrow(x)
  }

  #generate reduced matrix
  simiMat.ori <- simiMat
  simiMat[upper.tri(simiMat)] <- 0
  diag(simiMat) <- 0
  x.re <- x[,!apply(simiMat, 2, function(n) any(abs(n) >= cutoff))]

  if (record) {
    #a matrix keeping track of the removed features
    mapReduce <- simiMat.ori
    diag(mapReduce) <- 0
    mapList <- lapply(colnames(x.re), function(i) colnames(mapReduce)[mapReduce[i,]>=cutoff])
    names(mapList) <- colnames(x.re)
  } else mapList = NULL

  return(list(reduced = x.re,
              mapReduce = mapList))
}

# Run glmnet with gaussian familiy outcome
runGlm <- function(X, y, method = "lasso", repeats=20, folds = 3, testRatio = NULL, lambda = "lambda.1se") {
  modelList <- list()
  lambdaList <- c()
  r2Train <- c()
  r2Test <- c()
  coefMat  <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)

  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }

  for (i in seq(repeats)) {
    if (!is.null(testRatio)) {
      testIdx <- sample(seq_along(y), length(y)*testRatio)
      X.test <- X[testIdx,]
      X.train <- X[-testIdx,]
      y.test <- y[testIdx]
      y.train <- y[-testIdx]
    } else {
      X.train <- X
      y.train <- y
    }

    vecFold <- mltools::folds(y.train,nfolds = folds, stratified = TRUE)

    #train model
    res <- cv.glmnet(X.train,y.train, type.measure = "mse",
                     foldid = vecFold, alpha = alpha, standardize = FALSE,
                     intercept = TRUE, family = "gaussian")
    lambdaList <- c(lambdaList, res[[lambda]])

    #calculate variance explained for training model
    y.pred <- predict(res, s = lambda, newx = X.train)[,1]
    varExp <- cor(y.train,y.pred)^2
    r2Train <- c(r2Train, varExp)

    modelList[[i]] <- res
    coefMat[,i] <- coef(res, s = lambda)[-1]

    #test model if testRatio is speficied
    if(!is.null(testRatio)) {
      y.pred <- predict(res, s = lambda, newx = X.test)
      varExp <- cor(as.vector(y.test),as.vector(y.pred))^2
      r2Test <- c(r2Test, varExp)
    }
  }
  list(modelList = modelList, lambdaList = lambdaList, r2Train = r2Train, coefMat = coefMat,
       r2Test = r2Test)
}


#get number of selected features from a LASSO model
nFeature <- function(lassoRes) {
  coefMat <- lassoRes$coefMat
  return(colSums(coefMat != 0))
}



runCamera <- function(exprMat, design, gmtFile, id = NULL,
                      contrast = ncol(design),  method = "camera", pCut = 0.05,
                      ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE) {

  #prepare indices
  if (is.null(id)) id <- rownames(exprMat)

  if (is.character(gmtFile)) {
    idx <- limma::ids2indices(piano::loadGSC(gmtFile)$gsc, id)
  } else {
    idx <- limma::ids2indices(gmtFile,id)
  }

  #run camera for fry
  if (method == "camera") {
    res <- limma::camera(exprMat, idx, design, contrast)
  } else if (method == "fry") {
    res <- limma::fry(exprMat, idx, design, contrast)
  }

  #plot enrichment results as bar plot

  plotTab <- res %>% rownames_to_column("Name")
  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
  plotTab <- plotTab %>%
    mutate(Direction= factor(Direction, levels =c("Down","Up"))) %>%
    arrange(desc(Direction),desc(PValue)) %>%
    mutate(Name = factor(Name, levels = Name))

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = res, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue), fill = Direction)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5) +
      scale_fill_manual(values = c(Up = colList[1], Down = colList[2])) +
      coord_flip() + xlab("") +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_full +
      theme(axis.text = element_text(size= 12))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.1))
    } else {
      p <- p + theme(legend.position = "right")
    }

    return(list(enrichTab = res, enrichPlot = p))
  }
}

#function to run fisher test for enrichment analysis
runFisher <- function (genes, reference, gmtFile, adj = "BH", verbose = FALSE, barCol = NULL, setName = "",
                       pCut = 0.05,ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE) {

  genesets <- piano::loadGSC(gmtFile)$gsc
  tab = lapply(1:length(genesets), function(i) {
    if (verbose == TRUE) {
      cat("processing term", i, names(genesets)[i], "\n")
    }
    reference = reference[!reference %in% genes]
    RinSet = sum(reference %in% genesets[[i]])
    RninSet = length(reference) - RinSet
    GinSet = sum(genes %in% genesets[[i]])
    GninSet = length(genes) - GinSet
    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = F)
    colnames(fmat) = c("inSet", "ninSet")
    rownames(fmat) = c("genes", "reference")
    fish = fisher.test(fmat, alternative = "greater")
    pval = fish$p.value
    inSet = RinSet + GinSet
    res = c(GinSet, inSet, pval)
    res
  })
  rtab = do.call("rbind", tab)
  rtab = data.frame(as.vector(names(genesets)), rtab)
  rtab = rtab[order(rtab[, 4]), ]
  colnames(rtab) = c("TermID", "genes", "all", "pval")
  padj = p.adjust(rtab[, 4], method = "BH")
  tab.out = data.frame(rtab, padj)

  plotTab <- tab.out %>% dplyr::rename(Name = TermID, PValue = pval, FDR = padj)

  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))

  plotTab <- plotTab %>%
    arrange(desc(PValue)) %>%
    mutate(Name = sprintf("%s(%s)",Name,genes)) %>%
    mutate(Name = factor(Name, levels = Name))

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  if (is.null(barCol)) barCol <- colList[1]

  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = tab.out, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue))) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5, fill = barCol) +
      coord_flip() + xlab(setName) +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_full +
      theme(axis.text = element_text(size= 12))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.1))
    } else {
      p <- p + theme(legend.position = "right")
    }

    return(list(enrichTab = tab.out, enrichPlot = p))
  }
}



runFGSEA <- function(limmaRes, signatures, stat = "t", name = "symbol",minSize=15,
                     maxSize =1000,nperm = 10000, pcut = 0.1, ifFDR = TRUE, showFDR = FALSE) {
  geneRanks <- limmaRes[[stat]]
  names(geneRanks) <- limmaRes[[name]]
  fgseaRes <- fgsea(pathways = sigList,
                    stats = geneRanks,
                    minSize=minSize,
                    maxSize=maxSize,
                    nperm=nperm)
  if(ifFDR) {
    plotSets <- filter(fgseaRes, padj <= pcut)$pathway
  } else {
    plotSets <- filter(fgseaRes, pval <= pcut)$pathway
  }

  #plot gsea plots
  plots <- lapply(plotSets, function(setName) {
    pval <- formatNum(filter(fgseaRes, pathway == setName)$pval, digits = 1)
    FDR <- formatNum(filter(fgseaRes, pathway == setName)$padj, digits = 1)
    nes <- filter(fgseaRes, pathway == setName)$NES
    NES <- format(nes, digits = 1,nsmall = 1)
    #showValue <- ifelse(showFDR, sprintf("italic(P)~'value = %s\nFDR = %s\nNormalized ES=%s'", pval, FDR, NES),
    #                    sprintf("italic(P)~'value = %s\nNormalized ES=%s'", pval, NES))
    pp <- fgsea::plotEnrichment(signatures[[setName]],
                         geneRanks)
    pp <- pp + ggtitle(setName) +
      ylab("Enrichment score (ES)") + xlab("Rank") +
      theme(plot.title = element_text(face = "bold", size = 15),
            plot.margin = margin(2,1,2,1,"cm"),
            axis.text = element_text(size=15),
            axis.title = element_text(size=16)
            )
    if (nes >0) {
      showValue <- sprintf("atop('          '~italic(P)~'=%s','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = Inf, y = Inf,
                          label = showValue, parse= TRUE,
                          hjust=1, vjust =1.3, size = 5)
    } else {
      showValue <- sprintf("atop(italic(P)~'=%s            ','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = 0, y = -Inf,
                          label = showValue, parse=TRUE,
                          hjust=0, vjust =-0.5, size = 5)
    }
    pp
  })
  names(plots) <- plotSets
  return(list(table = fgseaRes, plots = plots))
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}

sumToTidy <- function(seObject, rowID = "rowID", colID = "colID") {

  tidyTable <- lapply(assayNames(seObject),function(n) {
    valTab <- assays(seObject)[[n]] %>% data.frame() %>%
      rownames_to_column(rowID) %>%
      gather(key = !!colID, value = "val", -!!rowID) %>%
      mutate(assay = n)
  }) %>% bind_rows() %>%
    spread(key = assay, value = val)

  #append row annotations
  rowAnno <- rowData(seObject) %>% data.frame() %>%
    rownames_to_column(rowID)
  tidyTable <- left_join(tidyTable, rowAnno, by = rowID)

  #append column annotations
  colAnno <- colData(seObject) %>% data.frame() %>%
    rownames_to_column(colID)
  tidyTable <- left_join(tidyTable, colAnno, by = colID)


  return(as_tibble(tidyTable))
}

runGSEA <- function(inputTab,gmtFile,GSAmethod="gsea",nPerm=1000){
  require(piano)
  inGMT <- loadGSC(gmtFile,type="gmt")
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] #re-rank by score

  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm, verbose = FALSE)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'nullDist', verbose = FALSE)
    GSAsummaryTable(res)
  }
}

plotEnrichmentBar <- function(resTab, pCut = 0.05, ifFDR = FALSE, setName = "", title="",
                              removePrefix = NULL, insideLegend = FALSE) {

    plotTab <- resTab

    if (ifFDR) {
      plotTab <- dplyr::filter(plotTab, `p adj (dist.dir.up)` <= pCut | `p adj (dist.dir.dn)` <= pCut)
    } else {
      plotTab <- dplyr::filter(plotTab, `p (dist.dir.up)` <= pCut | `p (dist.dir.dn)` <= pCut)
    }

    if (nrow(plotTab) == 0) {
      print("No sets passed the criteria")
      return(NULL)

    } else {
      #firstly, process the result table
      plotTab <- lapply(seq(nrow(plotTab)), function(i) {
        x <- plotTab[i,]
        statSign <- as.numeric(x[3])
        data.frame(Name = x[[1]], p = as.numeric(ifelse(statSign >= 0, x[[4]], x[[6]])),
                   geneNum = ifelse(statSign >= 0, x[[8]], x[[9]]),
                   Direction = ifelse(statSign > 0, "Up", "Down"), stringsAsFactors = FALSE)
      }) %>% bind_rows()

      if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))

      plotTab$Name <- sprintf("%s (%s)",plotTab$Name,plotTab$geneNum)
      plotTab <- plotTab[with(plotTab,order(Direction, p, decreasing=TRUE)),]
      plotTab$Direction <- factor(plotTab$Direction, levels = c("Up","Down"))
      plotTab$Name <- factor(plotTab$Name, levels = plotTab$Name)
      #plot the barplot
      p <- ggplot(data=plotTab, aes(x=Name, y= -log10(p), fill=Direction)) +
        geom_bar(position="dodge",stat="identity", width = 0.5) +
        scale_fill_manual(values=c(Up = colList[1], Down = colList[2])) +
        coord_flip() + xlab(setName) +
        ylab(expression(-log[10]*'('*p*')')) +
        ggtitle(title) + theme_full + theme(plot.title = element_text(face = "bold", hjust =0.5),
                                        axis.title = element_text(size=15))

      if (insideLegend) {
        p <- p + theme(legend.position = c(0.8,0.1))
      } else {
        p <- p + theme(legend.position = "right")
      }
    }


  return(p)
}



plotVolcano <- function(pTab, fdrCut = 0.05, posCol = "red", negCol = "blue",
                        x_lab = "dm", plotTitle = "",ifLabel = FALSE, labelList=NULL,
                        colLabel = NULL) {
  plotTab <- pTab %>% mutate(ifSig = ifelse(adj.P.Val > fdrCut, "n.s.",
                                            ifelse(logFC > 0, "up","down"))) %>%
    mutate(ifSig = factor(ifSig, levels = c("up","down","n.s.")))
  pCut <- -log10((filter(plotTab, ifSig != "n.s.") %>% arrange(desc(P.Value)))$P.Value[1])
  g <- ggplot(plotTab, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(shape = 21, aes(fill = ifSig),size=3) +
    geom_hline(yintercept = pCut, linetype = "dashed") +
    annotate("text", x = -Inf, y = pCut, label = paste0(fdrCut*100,"% FDR"),
             size = 5, vjust = -1.2, hjust=-0.1) +
    scale_fill_manual(values = c(n.s. = "grey70",
                                 up = posCol, down = negCol),name = "") +
    theme_full +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15)) +
    ylab(expression(-log[10]*'('*italic(P)~value*')')) +
    xlab(x_lab) + ggtitle(plotTitle)

  if (ifLabel) {
    if (is.null(labelList)) {
      labelTab <- plotTab
    } else {
      labelTab <- filter(plotTab, name %in% labelList, ifSig != "n.s.")
    }

    if (is.null(colLabel)) {
      g <- g + ggrepel::geom_label_repel(data = labelTab, aes(label = name),
                                         size=8, force = 5, col = "red")
    } else {
      g <- g+ggrepel::geom_label_repel(data = labelTab,
                                       aes_string(label = "name", col = colLabel),
                                       size=8, force = 5) +
        scale_color_manual(values = c(yes = "red",no = "black"))
    }
  }

  return(g)
}

plotBox <- function(plotTab, pValTabel = NULL, y_lab = "Protein expression") {
  plotList <- lapply(unique(plotTab$name), function(n) {
    eachTab <- filter(plotTab, !is.na(status), !is.na(count), name == n) %>%
      group_by(status) %>% mutate(n=n()) %>% ungroup() %>%
      mutate(group = sprintf("%s\n(N=%s)",status,n)) %>%
      arrange(status) %>% mutate(group = factor(group, levels = unique(group)))

    if (!is.null(pValTabel)) {
      pval <- formatNum(filter(pValTabel, name == n)$P.Value, digits = 1, format="e")
      titleText <- bquote(.(n)~" ("~italic("P")~"="~.(pval)~")")
    } else {
      titleText <- n
    }
    ggplot(eachTab, aes(x=group, y = count)) +
      geom_boxplot(width=0.3, aes(fill = group), outlier.shape = NA) +
      geom_beeswarm(col = "black", size =1,cex = 1.5, alpha=0.5) +
      ggtitle(titleText)+
      #ggtitle(sprintf("%s (p = %s)",geneName, formatNum(pval, digits = 1, format = "e"))) +
      ylab(y_lab) + xlab("") +
      scale_fill_manual(values = colList[3:5]) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
      theme_full +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.margin = margin(10,10,10,10))

  })
  return(plotList)
}

plotPep <- function(pepCLL, protName, type = "count", stratifier = NULL, mutStatus = NULL, protExpr = NULL) {

  #get peptide expression
  pepTab <- assays(pepCLL[rowData(pepCLL)$symbol %in% protName,])[[type]] %>%
    data.frame() %>% rownames_to_column("pepID") %>%
    gather(key = "patID", value = "expr", -pepID) %>%
    mutate(sequence = rowData(pepCLL[pepID,])$pepSeq) %>%
    dplyr::select(patID, expr, sequence) %>%
    mutate(measure = "Peptide")

  #prepare protein expression annotation table
  if (!is.null(protExpr)) {
    protTab <- tibble(patID = names(protExpr), expr = protExpr, sequence = "Protein", measure = "Protein")
    pepTab <- bind_rows(pepTab, protTab)
  }

  # prepare plotting table
  plotTab <- pepTab %>%
    group_by(patID) %>% mutate(medVal = median(expr, na.rm=TRUE)) %>%
    ungroup() %>% arrange(desc(medVal)) %>%
    mutate(patID = factor(patID, levels = unique(patID))) %>%
    dplyr::select(-medVal)

  #add genomic
  if (!is.null(mutStatus)) {
    plotTab <- mutate(plotTab, mut = mutStatus[match(patID, names(mutStatus))]) %>%
      filter(!is.na(mut)) %>% mutate(mut = paste0(stratifier, " = ", mut))
  }

  p <- ggplot(plotTab, aes(x = patID, y = expr, col = sequence, group = sequence)) +
    geom_point() +  theme_full +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab("") + ylab("Expression") + ggtitle(protName)

  #if plot protein expression level
  if (!is.null(protExpr)) {
    p <- p + geom_line(aes(linetype = measure)) + scale_shape_manual(values = c(Peptide = "solid",Protein = "dotted"))
  } else {
    p <- p + geom_line()
  }

  if (!is.null(mutStatus)) {
    p <- p + facet_wrap(~mut, scales = "free_x",ncol =length(unique(mutStatus)), drop = TRUE)
    gp <- ggplotGrob(p)
    facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
    x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                    function(l) length(l$range$range))
    gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
    p <- gp
  }

  p
}

compareR2 <- function(Drug, feature, geneMat, viabMat, dataMat) {
  viab <- viabMat[Drug,]
  expr <- dataMat[feature,]

  tabGene <- data.frame(geneMat)
  tabGene[["viab"]] <- viab
  tabCom <- tabGene
  tabCom[[feature]] <- expr

  r2Prot <- summary(lm(viab~expr))$adj.r.squared
  r2Gene <- summary(lm(viab~., data=tabGene))$adj.r.squared
  resCom <- summary(lm(viab~., data=tabCom))
  r2Com <- resCom$adj.r.squared
  pProt <- resCom$coefficients[feature,"Pr(>|t|)"]
  resTab <- tibble(Drug = Drug, feature = feature, r2Feature = r2Prot, r2Gene = r2Gene, r2Com = r2Com, r2Diff = r2Com - r2Gene,
                   p.value = pProt)
  resTab
}



runGlm <- function(X, y, method = "ridge", repeats=20, folds = 3) {
  modelList <- list()
  lambdaList <- c()
  varExplain <- c()
  coefMat <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)

  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }

  for (i in seq(repeats)) {
    if (ncol(X) > 2) {
      res <- cv.glmnet(X,y, type.measure = "mse", family="gaussian",
                       nfolds = folds, alpha = alpha, standardize = FALSE)
      lambdaList <- c(lambdaList, res$lambda.min)
      modelList[[i]] <- res

      coefModel <- coef(res, s = "lambda.min")[-1] #remove intercept row
      coefMat[,i] <- coefModel

      #calculate variance explained
      y.pred <- predict(res, s = "lambda.min", newx = X)
      varExp <- cor(as.vector(y),as.vector(y.pred))^2
      varExplain[i] <- ifelse(is.na(varExp), 0, varExp)

    } else {
      fitlm<-lm(y~., data.frame(X))
      varExp <- summary(fitlm)$r.squared
      varExplain <- c(varExplain, varExp)

    }

  }
  list(modelList = modelList, lambdaList = lambdaList, varExplain = varExplain, coefMat = coefMat)
}

plotVar <- function(glmResult, maxY =1) {

  pList <- lapply(names(glmResult), function(n) {
    plotTab <- lapply(names(glmResult[[n]]), function(x) {
      tibble(variable = x, value = glmResult[[n]][[x]]$varExplain)}) %>%
      bind_rows() %>% group_by(variable) %>%
      summarise(mean=mean(value, na.rm = TRUE),sd=sd(value, na.rm=TRUE))

    ggplot(plotTab,aes(x=variable, y=mean, fill= variable)) +
      geom_bar(position=position_dodge(), stat="identity", width = 0.8, col="black") +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, width = 0.3), position=position_dodge(.9)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =16),
                              plot.title = element_text(size=16, face ="bold"),
                              axis.text.y = element_text(size=16),
                              axis.title = element_text(size=16),
                              legend.position = "none") +
      scale_fill_brewer("Set1",type = "qual") + coord_cartesian(ylim = c(0,maxY)) +
      ylab("R2") + xlab("") + ggtitle(paste0("Variance explained for\n",n))
  })
  pList
}

lassoPlot <- function(lassoOut, cleanData, freqCut = 1, coefCut = 0.01, setNumber = "last") {
  library(gtable)
  plotList <- list()
  if (setNumber == "last") {
    setNumber <- length(lassoOut[[1]])
  } else {
    setNumber <- setNumber
  }
  for (seaName in names(lassoOut)) {
    #for the barplot on the left of the heatmap
    barValue <- rowMeans(lassoOut[[seaName]][[setNumber]]$coefMat)
    freqValue <- rowMeans(abs(sign(lassoOut[[seaName]][[setNumber]]$coefMat)))
    barValue <- barValue[abs(barValue) >= coefCut & freqValue >= freqCut] # a certain threshold
    barValue <- barValue[order(barValue)]
    if(length(barValue) == 0) {
      plotList[[seaName]] <- NA
      next
    }

    #for the heatmap and scatter plot below the heatmap
    allData <- cleanData$allExplain[[seaName]][[setNumber]]
    seaValue <- cleanData$allResponse[[seaName]]*2 #back to Z-score

    tabValue <- allData[, names(barValue),drop=FALSE]
    ord <- order(seaValue)
    seaValue <- seaValue[ord]
    tabValue <- tabValue[ord, ,drop=FALSE]
    sampleIDs <- rownames(tabValue)
    tabValue <- as.tibble(tabValue)

    #change scaled binary back to catagorical
    for (eachCol in colnames(tabValue)) {
      if (strsplit(eachCol, split = "[.]")[[1]][1] != "con") {
        tabValue[[eachCol]] <- as.integer(as.factor(tabValue[[eachCol]]))
      }
      else {
        tabValue[[eachCol]] <- tabValue[[eachCol]]*2 #back to Z-score
      }
    }

    tabValue$Sample <- sampleIDs
    #Mark different rows for different scaling in heatmap
    matValue <- gather(tabValue, key = "Var",value = "Value", -Sample)
    matValue$Type <- "mut"

    #For continuious value
    matValue$Type[grep("con.",matValue$Var)] <- "con"

    #for methylation_cluster
    matValue$Type[grep("ConsCluster",matValue$Var)] <- "meth"

    #change the scale of the value, let them do not overlap with each other
    matValue[matValue$Type == "mut",]$Value = matValue[matValue$Type == "mut",]$Value + 10
    matValue[matValue$Type == "meth",]$Value = matValue[matValue$Type == "meth",]$Value + 20


    #color scale for viability
    idx <- matValue$Type == "con"

    myCol <- colorRampPalette(c('dark red','white','dark blue'),
                              space = "Lab")
    if (sum(idx) != 0) {
      matValue[idx,]$Value = round(matValue[idx,]$Value,digits = 2)
      minViab <- min(matValue[idx,]$Value)
      maxViab <- max(matValue[idx,]$Value)
      limViab <- max(c(abs(minViab), abs(maxViab)))
      scaleSeq1 <- round(seq(-limViab, limViab,0.01), digits=2)
      color4viab <- setNames(myCol(length(scaleSeq1+1)), nm=scaleSeq1)
    } else {
      scaleSeq1 <- round(seq(0,1,0.01), digits=2)
      color4viab <- setNames(myCol(length(scaleSeq1+1)), nm=scaleSeq1)
    }



    #change continues measurement to discrete measurement
    matValue$Value <- factor(matValue$Value,levels = sort(unique(matValue$Value)))

    #change order of heatmap
    names(barValue) <-  gsub("con.", "", names(barValue))
    matValue$Var <- gsub("con.","",matValue$Var)
    matValue$Var <- factor(matValue$Var, levels = names(barValue))
    matValue$Sample <- factor(matValue$Sample, levels = names(seaValue))
    #plot the heatmap
    p1 <- ggplot(matValue, aes(x=Sample, y=Var)) + geom_tile(aes(fill=Value), color = "gray") +
      theme_bw() + scale_y_discrete(expand=c(0,0)) + theme(axis.text.y=element_text(hjust=0, size=12), axis.text.x=element_blank(), axis.ticks=element_blank(), panel.border=element_rect(colour="gainsboro"),
                                                           plot.title=element_text(size=14, face="bold"), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + xlab("patients") + ylab("") + scale_fill_manual(name="Mutated", values=c(color4viab, `11`="gray96", `12`='black', `21`='lightgreen', `22`='green',`23` = 'green4'),guide=FALSE) + ggtitle(seaName)

    #Plot the bar plot on the left of the heatmap
    barDF = data.frame(barValue, nm=factor(names(barValue),levels=names(barValue)))

    p2 <- ggplot(data=barDF, aes(x=nm, y=barValue)) +
      geom_bar(stat="identity", fill="lightblue", colour="black", position = "identity", width=.66, size=0.2) +
      theme_bw() + geom_hline(yintercept=0, size=0.3) + scale_x_discrete(expand=c(0,0.5)) +
      scale_y_continuous(expand=c(0,0), breaks = scales::breaks_extended(n = 4)) + coord_flip() +
      theme(panel.grid.major=element_blank(), panel.background=element_blank(), axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank(), axis.text=element_text(size=10), panel.border=element_blank()) +
      xlab("") + ylab("") + geom_vline(xintercept=c(0.5), color="black", size=0.6)

    #Plot the scatter plot under the heatmap

    # scatterplot below
    scatterDF = data.frame(X=factor(names(seaValue), levels=names(seaValue)), Y=seaValue)

    p3 <- ggplot(scatterDF, aes(x=X, y=Y)) + geom_point(shape=21, fill="dimgrey", colour="black", size=1.2) + theme_bw() +
      theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), axis.text.y=element_text(size=10), panel.border=element_rect(colour="dimgrey", size=0.1), panel.background=element_rect(fill="gray96"))

    #Scale bar for continuous variable

    Vgg = ggplot(data=data.frame(x=1, y=as.numeric(names(color4viab))), aes(x=x, y=y, color=y)) + geom_point() +
      scale_color_gradientn(name="Z-score", colours =color4viab) + theme(legend.title=element_text(size=12), legend.text=element_text(size=10))

    #Assemble all the plots togehter

    # construct the gtable
    wdths = c(1.5, 0.2, 1.1*ncol(matValue), 0.9)
    hghts = c(0.1, 0.0025*nrow(matValue), 0.12, 0.6)
    gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))

    ## make grobs
    gg1 = ggplotGrob(p1)
    gg2 = ggplotGrob(p2)
    gg3 = ggplotGrob(p3)
    #gg4 = ggplotGrob(Vgg)

    ## fill in the gtable
    gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 2, 1) # add barplot
    gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 2, 3) # add heatmap
    gt = gtable_add_grob(gt, gtable_filter(gg1, "title"), 1, 3) #add title to plot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "panel"), 4, 3) # add scatterplot
    gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-b"), 3, 1) # y axis for barplot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "axis-l"), 4, 2) # y axis for scatter plot
    gt = gtable_add_grob(gt, gtable_filter(gg1, "axis-l"), 2, 4) # variable names
    #gt = gtable_add_grob(gt, gtable_filter(gg4, "guide-box"), 2, 5) # scale bar for continous variables


    #plot
    #grid.draw(gt)
    plotList[[seaName]] <- gt
  }
  return(plotList)
}

dataScale <- function(x, censor = NULL, robust = FALSE) {
  #function to scale different variables
  if (length(unique(na.omit(x))) == 2){
    #a binary variable, change to -0.5 and 0.5 for 1 and 2
    x - 0.5
  } else if (length(unique(na.omit(x))) == 3) {
    #catagorical varialbe with 3 levels, methylation_cluster, change to -0.5,0,0.5
    (x - 1)/2
  } else {
    if (robust) {
      #continuous variable, centered by median and divied by 2*mad
      mScore <- (x-median(x,na.rm=TRUE))/mad(x,na.rm=TRUE)
      if (!is.null(censor)) {
        mScore[mScore > censor] <- censor
        mScore[mScore < -censor] <- -censor
      }
      mScore/2
    } else {
      mScore <- (x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE))
      if (!is.null(censor)) {
        mScore[mScore > censor] <- censor
        mScore[mScore < -censor] <- -censor
      }
      mScore/2
    }
  }
}

#function to generate response vector and explainatory variable for each seahorse measurement
generateData.drug <- function(inclSet, onlyCombine = FALSE, censor = NULL, robust = FALSE) {


  allResponse <- list()
  allExplain <- list()

  for (measure in colnames(inclSet$drug)) {
    y <- inclSet$drug[,measure]
    y <- y[!is.na(y)]

    #get overlapped samples for each dataset
    overSample <- names(y)

    for (eachSet in inclSet) {
      overSample <- intersect(overSample,rownames(eachSet))
    }

    y <- dataScale(y[overSample], censor = censor, robust = robust)

    #generate explainatory variable table for each seahorse measurement
    expTab <- list()


    if ("gen" %in% names(inclSet)) {
      geneTab <- inclSet$gen[overSample,]
      #at least 3 mutated sample
      geneTab <- geneTab[, colSums(geneTab) >= 3]
      vecName <- sprintf("genetic(%s)", ncol(geneTab))
      expTab[[vecName]] <- apply(geneTab,2,dataScale)
    }

    if ("BH3" %in% names(inclSet)){

      #for BH3 profling
      bh3Data <- inclSet$BH3[overSample, ]
      colnames(bh3Data) <- paste0("con.",colnames(bh3Data), sep = "")
      vecName <- sprintf("BH3(%s)", ncol(bh3Data))
      expTab[[vecName]] <- apply(bh3Data,2,dataScale, censor = censor, robust = robust)

    }

    comboTab <- c()
    for (eachSet in names(expTab)){
      comboTab <- cbind(comboTab, expTab[[eachSet]])
    }
    vecName <- sprintf("all(%s)", ncol(comboTab))
    expTab[[vecName]] <- comboTab

    allResponse[[measure]] <- y
    allExplain[[measure]] <- expTab
  }
  if (onlyCombine) {
    #only return combined results, for feature selection
    allExplain <- lapply(allExplain, function(x) x[length(x)])
  }

  return(list(allResponse=allResponse, allExplain=allExplain))

}
