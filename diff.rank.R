diff.rank <- function(reverse.pair = NULL, expr = NULL, tum.sam = NULL, nor.sam = NULL) {
  
  expr <- expr[unique(unlist(strsplit(reverse.pair,"-"))),]
  
  expr.rank <- apply(expr, 2, rank)
  
  diff.rank.outTab <- NULL
  for (i in reverse.pair) {
    g1 <- unlist(strsplit(i,"-"))[1]
    g2 <- unlist(strsplit(i,"-"))[2]
    
    tmp <- expr.rank[g1,] - expr.rank[g2,]
    tmp <- as.data.frame(t(tmp)); rownames(tmp) <- i
    diff.rank.outTab <- rbind.data.frame(diff.rank.outTab,tmp)
  }
  
  diff.rank.tum <- diff.rank.outTab[,tum.sam]
  diff.rank.nor <- diff.rank.outTab[,nor.sam]
  
  diff.rank.tum <- as.data.frame(abs(pmax(as.matrix(diff.rank.tum),0)))
  diff.rank.nor <- as.data.frame(abs(pmin(as.matrix(diff.rank.nor),0)))
  
  return(list(diff.rank.tum = diff.rank.tum, diff.rank.nor = diff.rank.nor))
}