slec_stable_pair_nor <- function(expE,cutoff) { 
  
  library(matlab)
  gid <- rownames(expE)
  exp <- as.matrix(expE)
  len <- length(gid)-1 
  end <- length(gid)
  m   <- length(exp[1,]) 
  pair <- matrix(nrow=choose(length(gid),2),ncol=3) 
  j <- 1 
  for (i in 1:len) {
    RE=(exp[i*ones(len-i+1,1),]-exp[(i+1):end,]) < 0 
    if(i == len) {
      ratio <- sum(RE)/m
      pair[j,] <- c(rep(gid[i],end-i),gid[(i+1):end],ratio) 
    } else {
      ratio=rowSums(RE)/m 
      pair[j:(j+length(RE[,1])-1),] <- c(rep(gid[i],end-i),gid[(i+1):end],ratio) 
      j <- j + length(RE[,1])
    }
    
  } 
  stable_pair_1 <- pair[pair[,3] >= cutoff,1:2] 
  stable_pair_2 <- pair[pair[,3] <= 1-cutoff,2:1] 
  stable_pair <- rbind.data.frame(stable_pair_1,stable_pair_2) 
  rownames(stable_pair) <- paste0(stable_pair$V1,"-",stable_pair$V2)
  return(stable_pair) 
} 