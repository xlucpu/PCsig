valuate.performance <- function(optimal.sig.pair = NULL, expr = NULL, TrueClass = NULL) {
  
  genes <- unique(unlist(strsplit(optimal.sig.pair,"-")))
  loss.genes <- setdiff(genes,rownames(expr))
  if(length(loss.genes) > 0) {
    optimal.sig.pair2 <- optimal.sig.pair[-grep(paste(loss.genes,collapse="|"),optimal.sig.pair)]
  } else {optimal.sig.pair2 <- optimal.sig.pair}
  
  for (p in optimal.sig.pair2) {
    
    g1 <- unlist(strsplit(p,"-"))[1]
    g2 <- unlist(strsplit(p,"-"))[2]
    
    if (p == optimal.sig.pair2[1]) {
      diff.expr <- data.frame(diff = apply(expr,2, function(x) x[g1]-x[g2]),stringsAsFactors = F); colnames(diff.expr) <- p
    } else {
      tmp <- data.frame(diff = apply(expr,2, function(x) x[g1]-x[g2]),stringsAsFactors = F); colnames(tmp) <- p
      diff.expr <- cbind.data.frame(diff.expr,tmp)
    }
  }
  diff.expr <- as.data.frame(ifelse(diff.expr > 0, 1, 0))
  diff.expr$count <- as.numeric(rowSums(diff.expr))
  # rownames(diff.expr)[which(is.na(diff.expr$count))]
  # all(diff.expr$count >=0)
  diff.expr$TrueClass <- TrueClass
  
  # ROC
  library(tdROC)
  TP.value <- c()
  FP.value <- c()
  FN.value <- c()
  TN.value <- c()
  TPR.value <- c()
  FPR.value <- c()
  
  #for (pb in seq(1,max(diff.expr$count),1)) {
  for (pb in seq(1,length(optimal.sig.pair2),1)) {
    
    tmp <- diff.expr
    tmp$PredClass <- ifelse(tmp$count >= pb,"Tumor","Normal")
    #if(pb == ifelse(max(diff.expr$count) == 1,1,ceiling(0.5*max(diff.expr$count)))) {
    cutoff <- ifelse(length(optimal.sig.pair2) %% 2 == 0,
                     ceiling(0.5*length(optimal.sig.pair2)) + 1,
                     ceiling(0.5*length(optimal.sig.pair2)))
    if(pb == cutoff) {
      general.cutoff <- tmp
    }
    TP = sum(tmp$PredClass == "Tumor" & tmp$TrueClass == "Tumor")
    FP = sum(tmp$PredClass == "Tumor" & tmp$TrueClass == "Normal")
    FN = sum(tmp$PredClass == "Normal" & tmp$TrueClass == "Tumor")
    TN = sum(tmp$PredClass == "Normal" & tmp$TrueClass == "Normal")
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    
    TP.value <- c(TP.value,TP)
    FP.value <- c(FP.value,FP)
    FN.value <- c(FN.value,FN)
    TN.value <- c(TN.value,TN)
    TPR.value <- c(TPR.value,TPR)
    FPR.value <- c(FPR.value,FPR)
  }
  ROC.point <- data.frame(FPR = FPR.value, TPR = TPR.value, Spec = 1-FPR.value, Sens = TPR.value)
  auc <- calc.AUC(1-ROC.point$FPR,ROC.point$TPR)
  
  return(list(diff.expr = diff.expr,general.cutoff=general.cutoff,loss.genes=loss.genes, ROC.point = ROC.point, AUC=auc))
}
