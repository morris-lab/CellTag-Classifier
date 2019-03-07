#############################################################################################################
# Name: CellTag Classification
# Description: This file contains the source code for carrying out celltag Classification
#
# Functions:
# - normalize.function
# - permutation.sampling
# - dynamic.celltag.detection
# - binarization.classification
# - multiplet.calling
# - multiplet.checkpoint
#
# Additional detailed information can be found in the following reference:
# - Guo et. al., CellTag Indexing: a genetic barcode-based multiplexing tool for single-cell technologies, 
#   bioRxiv, 2018
#############################################################################################################

#############################################################################################################
# Function: normalize.function
# Description: This function does normalizations for CellTag matrices.
#
# Input
# - raw.log.mtx: logged matrix from count matrices
# - ignore.all.zero.column: would you like to ignore the columns with all zeros for size factor calculation
# Output
# - ct.dge.norm: a normalized matrix
#############################################################################################################
normalize.function <- function(raw.log.mtx, ignore.all.zero.column = TRUE) {
  ct.dge.norm <- raw.log.mtx
  csums <- colSums(raw.log.mtx)
  csums.avg <- mean(csums)
  if (ignore.all.zero.column) {
    csums.avg <- mean(csums[which(csums > 0)])
  }
  for (k in 1:length(csums)) {
    if (csums[k] == 0) {
      ct.dge.norm[,k] <- 0
    } else {
      ct.dge.norm[,k] <- (ct.dge.norm[,k]/csums[k]) * csums.avg
    }
  }
  return(ct.dge.norm)
}

#############################################################################################################
# Function: permutation.sampling
# Description: This function randomly draw given number of samples from the whole dataset and calculate
#              the likelihood of getting the benchmarking value or higher.
#
# Input
# - x: source dataset to draw samples from
# - curr.val: benchmarking value
# - n: sample size. Default is 1000 samples.
# Output
# - perc: a number that represents the likelihood of getting the benchmarking value or higher
#############################################################################################################
permutation.sampling <- function(x, curr.val, n = 1000) {
  samp <- sample(x, n, replace = T)
  perc <- sum(samp >= curr.val)/n
  return(perc)
}

#############################################################################################################
# Function: dynamic.celltag.detection
# Description: This function dynamically calculate the likelihood of occurrence of each celltag expression for
#              given number of iterations.
#
# Input
# - ct.dge.norm: normalized CellTag expression matrix. Expecting the matrix having row = cell, column = gene.
# - sample.num: sample size to draw. Default is 1000.
# - test.times: number of iterations to draw the samples. Default is 50.
# Output
# - perc.ls: a list contains the permutation probability results for each CellTag in each Cell.
#
# Note:
# 1. Each term of the output list contains the following information
#    a. Cell barcode
#    b. Matrix containing permutation probabilities
# 2. This function could take a long time!
#############################################################################################################
dynamic.celltag.detection <- function(ct.dge.norm, sample.num = 1000, test.times = 50){
  # Set up the initial values
  perc.ls <- list()
  count <- 1
  
  # Calculate the permutation probabilities
  for (l in 1:nrow(ct.dge.norm)) {
    curr.df <- data.frame()
    curr.cell <- rownames(ct.dge.norm)[l]
    for (m in 1:ncol(ct.dge.norm)){
      curr.col <- ct.dge.norm[,m]
      dens <- density(curr.col, n = length(curr.col))
      
      curr.num <- ct.dge.norm[l,m]
      perc.vec <- c()
      if (curr.num == 0) {
        perc.vec <- rep(1, test.times)
      } else {
        for (h in 1:test.times) {
          perc.vec[h] <- permutation.sampling(dens$x, curr.num, n = sample.num)
        }
      }
      if (nrow(curr.df) <= 0) {
        curr.df <- data.frame(perc.vec)
      } else {
        curr.df <- cbind(curr.df, perc.vec)
      }
      colnames(curr.df)[m] <- colnames(ct.dge.norm)[m]
    }
    perc.ls[[count]] <- list(curr.cell, curr.df)
    count <- count + 1
  }
  return(perc.ls)
}

#############################################################################################################
# Function: binarization.classification
# Description: This function uses the percentage list calculated from dynamic celltag detection to 
#              binarize and classify CellTag existences in each cell to 0 and 1.
#
# Input
# - ct.dge.norm: normalized CellTag expression matrix. Expecting the matrix having row = cell, column = gene.
# - perc.list: permutation probability list for all cells
# - multiplet.diff: differences used during pair-wise comparison of tags to call multiplet
# Output
# - bin.class.ct: a binary matrix of each cell vs. each celltag
#
# Note:
# 1. binary matrix:
#    - 0: should not be classified as this tag
#    - 1: should be classified as this tag
#############################################################################################################
binarization.classification <- function(ct.dge.norm, perc.list, multiplet.diff = 0.238){
  bin.class.ct <- data.frame(row.names = rownames(ct.dge.norm), stringsAsFactors = F)
  bin.class.ct[,colnames(ct.dge.norm)] <- 0
  
  for (i in 1:length(perc.list)) {
    curr.cell <- perc.list[[i]][[1]]
    curr.perc.df <- perc.list[[i]][[2]]
    cmeans <- colMeans(curr.perc.df)
    if (sum(cmeans == 1) == ncol(curr.perc.df)) {
      bin.class.ct[curr.cell, ] <- 0
    } else {
      min.col.num <- apply(perc.list[[i]][[2]], 1, function(x) which(x == min(x))[1])
      col.n.uniq <- unique(min.col.num)
      if (length(col.n.uniq) == 1) {
        min.col <- curr.perc.df[, col.n.uniq]
        other.perc.col <- as.data.frame(curr.perc.df[,-c(col.n.uniq, which(cmeans == 1))])
        colnames(other.perc.col) <- colnames(curr.perc.df)[-c(col.n.uniq, which(cmeans == 1))]
        other.perc.fold <- abs(min.col - other.perc.col)
        appx.one.rate <- apply(as.data.frame(other.perc.fold), 2, median)
        multi.col <- names(appx.one.rate[which(appx.one.rate <= multiplet.diff)])
        if (length(multi.col) >= 1) {
          bin.class.ct[curr.cell, multi.col] <- 1
        }
      }
      if (sum(cmeans[col.n.uniq] > 0.7) == length(col.n.uniq)) {bin.class.ct[curr.cell, ] <- 0}
      else {bin.class.ct[curr.cell, col.n.uniq] <- 1}
      bin.class.ct[curr.cell, which(cmeans == 1)] <- 0
    }
  }
  return(bin.class.ct)
}

#############################################################################################################
# Function: multiplet.calling
# Description: This function calls multiplets based on the binarized matrix
#
# Input
# - curr.ct.binary: binarized CellTag called matrix
# 
# Output
# - bin.class.ct: a binary & CellTag called matrix
#
# Note:
# 1. binary matrix:
#    - 0: should not be classified as this tag
#    - 1: should be classified as this tag
#    - last column: Multiplet/<Single Tag>
#############################################################################################################
multiplet.calling <- function(curr.ct.binary) {
  # assign classification based on binarized celltag expression
  bin.class.ct <- curr.ct.binary
  key <- colnames(curr.ct.binary)
  for (p in 1:nrow(bin.class.ct)) {
    if (sum(bin.class.ct[p,c(1:length(key))]) == 0) {
      bin.class.ct[p, "ct.call"] = "nd"
    } else if (sum(bin.class.ct[p,c(1:length(key))]) == 1) {
      bin.class.ct[p, "ct.call"] = key[which(bin.class.ct[p,] == 1)]
    } else if (sum(bin.class.ct[p,c(1:length(key))]) > 1) {
      bin.class.ct[p, "ct.call"] = "multiplet"
    }
  }
  return(bin.class.ct)
}

#############################################################################################################
# Function: multiplet.checkpoint
# Description: This function extracts the identified multiplet and recheck the multiplet to eliminate filter the 
#              parts that exceed the expected multiplet number
#
# Input
# - curr.ct.binary: binarized celltag called matrix
# - perc.list: permutation probability list for all cells
# - multiplet.table.path: where did you store your multiplet rate table?
# Output
# - bin.class.ct: a filtered binary CellTag called matrix
#
# Note:
# 1. binary matrix:
#    - 0: should not be classified as this tag
#    - 1: should be classified as this tag
#    - last column: Multiplet/<single tag>
#############################################################################################################
multiplet.checkpoint <- function(curr.ct.binary, perc.list, multiplet.table.path){
  bin.class.ct <- curr.ct.binary
  
  # Set up multiplet rate table
  multi.table <- read.csv(multiplet.table.path, header = T, stringsAsFactors = F)
  x <- multi.table$recovered.cell.num
  y <- multi.table$rate
  fit <- lm(y~x)
  stdev <- sqrt(var(fit$residuals))
  
  # Set up the expected multiplet rate for this dataset
  expected.multi.rate <- predict(fit, newdata = data.frame(x=nrow(bin.class.ct)))
  expected.multi.num <- nrow(bin.class.ct) * expected.multi.rate/100
  std.num <- nrow(bin.class.ct) * stdev/100
  
  multi.subset <- bin.class.ct[which(bin.class.ct$ct.call == "multiplet"), c(1:(ncol(bin.class.ct) - 1))]
  new.perc.inspec.ls <- list()
  if (nrow(multi.subset) > (expected.multi.num + std.num)) {
    for (i in 1:nrow(multi.subset)) {
      curr.cell <- rownames(multi.subset)[i]
      indx <- which(rownames(bin.class.ct) == curr.cell)
      new.perc.inspec.ls[[i]] <- perc.list[[indx]]
      new.perc.inspec.ls[[i]][[3]] <- colnames(multi.subset)[which(multi.subset[curr.cell,] == 1)]
    }
    min.diff.ls <- lapply(new.perc.inspec.ls, 
                          function(x) {
                            curr.col.subset <- x[[2]][, x[[3]]]
                            min.col.num <- apply(curr.col.subset, 1, function(x) which(x == min(x))[1])
                            col.n.uniq <- unique(min.col.num)
                            if (length(col.n.uniq) > 1) {
                              return(list(x[[1]], 0,0, NA))
                            } else {
                              curr.diff.mean.vec <- c()
                              for (i in 1:ncol(curr.col.subset)) {
                                if (i != col.n.uniq) {
                                  curr.diff <- mean(abs(curr.col.subset[,i] - curr.col.subset[,col.n.uniq]))
                                  curr.diff.mean.vec <- c(curr.diff.mean.vec, curr.diff)
                                }
                              }
                              return(list(x[[1]], min(curr.diff.mean.vec),max(curr.diff.mean.vec), x[[3]][col.n.uniq]))
                            }
                          })
    
    min.diff.df <- data.frame()
    for (j in 1:length(min.diff.ls)) {
      curr.list <- min.diff.ls[[j]]
      curr.df <- data.frame(cell.bc = curr.list[[1]], min.diff = curr.list[[2]],max.diff = curr.list[[3]], unique.col = curr.list[[4]], stringsAsFactors = F)
      if (nrow(min.diff.df) <= 0) {
        min.diff.df <- curr.df
      } else {
        min.diff.df <- rbind(min.diff.df, curr.df)
      }
    }
    min.diff.df.sort <- min.diff.df[order(-min.diff.df$min.diff), ]
    min.diff.df.abv <- min.diff.df[which(min.diff.df$min.diff > (quantile(min.diff.df$min.diff, 1.5 * expected.multi.num/nrow(min.diff.df.sort)))), ]
    
    for (k in 1:nrow(min.diff.df.abv)) {
      curr.cell.bc <- min.diff.df.abv$cell.bc[k]
      curr.tag.to.set <- min.diff.df.abv$unique.col[k]
      indx.to.erase <- which(colnames(multi.subset) != curr.tag.to.set)
      bin.class.ct[curr.cell.bc, indx.to.erase] <- 0
      bin.class.ct[curr.cell.bc, "ct.call"] <- curr.tag.to.set
    }
  } else {
    print(paste0("Multiplet amount not exceeding the expected number = ", round(expected.multi.num), ". Secondary check not applied!"))
  }
  return(bin.class.ct)
}





