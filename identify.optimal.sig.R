identify.optimal.sig <- function(reverse.pair.avgR = NULL, expr = NULL, TrueClass = NULL) {
  
  # library(tdROC)
  TP.value  <- c()
  FP.value  <- c()
  FN.value  <- c()
  TN.value  <- c()
  TPR.value <- c()
  FPR.value <- c()
  ACC.value <- c()
  #auc.value <- c()
  
  for (i in 1:nrow(reverse.pair.avgR)) {
    
    gene.pair <- reverse.pair.avgR$pair[1:i]
    
    diff.expr <- NULL
    for (p in gene.pair) {
      g1 <- unlist(strsplit(p,"-"))[1]
      g2 <- unlist(strsplit(p,"-"))[2]
      
      if (p == gene.pair[1]) {
        diff.expr <- data.frame(diff = apply(expr,2, function(x) x[g1]-x[g2]),stringsAsFactors = F); colnames(diff.expr) <- p
      } else {
        tmp <- data.frame(diff = apply(expr,2, function(x) x[g1]-x[g2]),stringsAsFactors = F); colnames(tmp) <- p
        diff.expr <- cbind.data.frame(diff.expr,tmp)
      }
    }
    diff.expr$PredClass <- ifelse(rowSums(diff.expr > 0) >= ceil(ncol(diff.expr)/2), "Tumor","Normal")
    diff.expr$TrueClass <- TrueClass
    
    TP = sum(diff.expr$PredClass == "Tumor" & diff.expr$TrueClass == "Tumor")
    FP = sum(diff.expr$PredClass == "Tumor" & diff.expr$TrueClass == "Normal")
    FN = sum(diff.expr$PredClass == "Normal" & diff.expr$TrueClass == "Tumor")
    TN = sum(diff.expr$PredClass == "Normal" & diff.expr$TrueClass == "Normal")
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    ACC = (TP + TN) / ncol(expr)
    
    TP.value <- c(TP.value,TP)
    FP.value <- c(FP.value,FP)
    FN.value <- c(FN.value,FN)
    TN.value <- c(TN.value,TN)
    TPR.value <- c(TPR.value,TPR)
    FPR.value <- c(FPR.value,FPR)
    ACC.value <- c(ACC.value,ACC)
  }
  
  # ROC.point <- data.frame(FPR = FPR.value, TPR = TPR.value)
  # auc.value <- calc.AUC(1-ROC.point$FPR,ROC.point$TPR)
  
  reverse.pair.avgR$Sensitivity <- TPR.value
  reverse.pair.avgR$Specificity <- 1-FPR.value
  reverse.pair.avgR$Accuracy <- ACC.value
  
  return(reverse.pair.avgR)
}
