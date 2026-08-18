###########################################################################################
#################### Create plots for manuscript (Dist to input) ##########################
###########################################################################################

### This script was written by a human and optimized using Claude LLM


library(openxlsx)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(MSCquartets)

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"
sim_out  <- file.path(base_dir, "sim_out")
source(paste0(base_dir, "/auxfiles/Finalcode_bootsupport.R"))

n_reps <- 10

measure_levels  <- c("Blobs","Cycles","Hybrids B>4","Splits","Mean")
category_rooted <- c("TRUE","Squirrel","HMS","HWC","L1N")
category_unroot <- c("TRUE","Squirrel","L1N-unrooted")
cols_rooted <- c("gray35","gray70","#A8DADC","#274C77","#A8D5A2")
cols_unroot <- c("gray35","gray70","#1B5E20")

 
big_text <- theme(
  plot.title   = element_text(size = 18, hjust = 0.5 ),
  axis.text    = element_text(size = 16),
  axis.title   = element_text(size = 18),
  legend.text  = element_text(size = 17),
  legend.title = element_text(size = 17),
  legend.key.size = unit(1.2, "lines")
)

 
parse_canon <- function(canon_str) {
  blocks <- strsplit(canon_str, " / ", fixed = TRUE)[[1]]
  lapply(blocks, function(b) strsplit(b, "*", fixed = TRUE)[[1]])
}

 
restrict_partition <- function(blocks, taxa_keep) {
  blocks <- lapply(blocks, function(b) sort(intersect(b, taxa_keep)))
  blocks <- blocks[lengths(blocks) > 0]
  blocks[order(vapply(blocks, `[`, character(1), 1))]
}

 
canon_compatible <- function(canon1, canon2) {
  b1 <- parse_canon(canon1)
  b2 <- parse_canon(canon2)
  common <- intersect(unlist(b1), unlist(b2))
  if (length(common) == 0) return(FALSE)
  identical(restrict_partition(b1, common), restrict_partition(b2, common))
}

 
merge_zero_rows <- function(df) {
  score_cols <- setdiff(names(df), c("canon", "n_partitions"))
  first_col <- score_cols[1]
  
  repeat {
    zero_rows <- which(df[[first_col]] == 0)
    if (length(zero_rows) == 0) break
    
    i <- zero_rows[which.min(df$n_partitions[zero_rows])]
    s <- df$n_partitions[i]
    
    nonzero_cols_i <- score_cols[as.numeric(df[i, score_cols]) != 0]
    if (length(nonzero_cols_i) == 0) break
    
    candidates <- which(
      df$n_partitions > s &
        apply(df[, nonzero_cols_i, drop = FALSE], 1,
              function(r) all(as.numeric(r) == 0)) &
        vapply(df$canon, canon_compatible, logical(1), canon2 = df$canon[i])
    )
    
    if (length(candidates) == 0) break
    
    j <- candidates[which.min(df$n_partitions[candidates])]
    
    for (col in nonzero_cols_i) {
      df[j, col] <- df[i, col]
    }
    df <- df[-i, ]
  }
  
  rownames(df) <- NULL
  return(df)
}
 
compute_NB_matrix <- function(location, netfile_fn) {
  NB <- matrix(rep(0, 2*20), nrow = 20)
  for (xx in 1:20) {
    if (xx < 10) { xp <- paste0('', xx) } else { xp <- xx }
    
    netfile <- file.path(location, netfile_fn(xp))
    sboot_sampled_lvl1 <- read.nets(file = netfile)[2:101]
    sboot_sampled_tob <- list()
    nblobs  <- 0
    n4blobs <- 0
    
    for (i in 1:100) {
      flag <- 1
      if (is.null(sboot_sampled_lvl1[[i]]$reticulation)) { flag <- 0 }
      if (flag == 1) {
        sboot_sampled_tob[[i]] <- treeOfBlobs(sboot_sampled_lvl1[[i]])
      } else {
        sboot_sampled_tob[[i]] <- sboot_sampled_lvl1[[i]]
      }
      deg <- ape::degree(sboot_sampled_tob[[i]])
      nblobs  <- sum(deg[deg[,1] > 3, 2]) + nblobs
      n4blobs <- sum(deg[deg[,1] > 4, 2]) + n4blobs
    }
    NB[xx, 1] <- nblobs
    NB[xx, 2] <- n4blobs
  }
  NB
}
 
compute_rooted_metrics <- function(NN, L, M) {
  
  t1blobs_a = rep(0,5); t1blobs_b = rep(0,5); t1blobs_c = rep(0,5)
  t1CO_a    = rep(0,5); t1CO_b    = rep(0,5); t1CO_c    = rep(0,5)
  t1Hy_a    = rep(0,5); t1Hy_b    = rep(0,5); t1Hy_c    = rep(0,5)
  t1Sp      = rep(0,5)
  t2blobs_a = rep(0,5); t2blobs_b = rep(0,5); t2blobs_c = rep(0,5)
  t2CO_a    = rep(0,5); t2CO_b    = rep(0,5); t2CO_c    = rep(0,5)
  t2Hy_a    = rep(0,5); t2Hy_b    = rep(0,5); t2Hy_c    = rep(0,5)
  t2Sp      = rep(0,5)
  
  t1Sp_a = rep(0,5); t1Sp_b = rep(0,5); t1Sp_c = rep(0,5)
  t2Sp_a = rep(0,5); t2Sp_b = rep(0,5); t2Sp_c = rep(0,5)
  
  t1Sp_a[1]=sum(M[,29]); t1Sp_b[1]=sum(1*M[,3]); t1Sp_c = rep(sum(1*M[,1]), 5)
  t1Sp_a[2]=sum(M[,28]); t1Sp_b[2]=sum(1*M[,2])
  t1Sp_a[3]=sum(M[,30]); t1Sp_b[3]=sum(1*M[,4])
  t1Sp_a[4]=sum(M[,31]); t1Sp_b[4]=sum(1*M[,5])
  t1Sp_a[5]=sum(M[,27]); t1Sp_b[5]=sum(1*M[,1])
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(L[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(L[[i]][[1]][L[[i]][[1]][, 2] > 4, 8])
    
    yy=0
    if(L[[i]][[yy+1]][1,1]=='NA'){
      flag=1
    } else {
      L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
      
      t1blobs_a[yy+5] = t1blobs_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][,4]))*100
      t1blobs_c[yy+5] = t1blobs_c[yy+5] + 1*truen_num_blobs*100
      t1blobs_b[yy+5] = t1blobs_b[yy+5] + 1*sum(L[[i]][[yy+1]][,3]>0)*100
      
      t2blobs_a[yy+5] = t2blobs_a[yy+5] +  sum(as.numeric(L[[i]][[yy+1]][,4]))*100
      t2blobs_c[yy+5] = t2blobs_c[yy+5] + 1*truen_num_blobs*100
      
      t1CO_a[yy+5] = t1CO_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][,6]))*100
      t1CO_c[yy+5] = t1CO_c[yy+5] + 1*truen_num_blobs*100
      t1CO_b[yy+5] = t1CO_b[yy+5] + 1*sum(L[[i]][[yy+1]][,6]>0)*100
      
      t2CO_a[yy+5] = t2CO_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][,6]))*100
      t2CO_c[yy+5] = t2CO_c[yy+5] + 1*truen_num_blobs*100
      
      if(truen_num_blobs_ge4 > 0)
      {
        t1Hy_a[yy+5] = t1Hy_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 8] ))*100
        t1Hy_c[yy+5] = t1Hy_c[yy+5] + 1*truen_num_blobs_ge4*100
        t1Hy_b[yy+5] = t1Hy_b[yy+5] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 7]) *100
        
        t2Hy_a[yy+5] = t2Hy_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 8] ))*100
        t2Hy_c[yy+5] = t2Hy_c[yy+5] + 1*truen_num_blobs_ge4*100
      }
    }
    
    for(yy in 1:2)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        flag=1
      } else {
        L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,4]))*100
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs*100
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0)*100
        
        t2blobs_a[yy] = t2blobs_a[yy] +  sum(as.numeric(L[[i]][[yy+1]][,4]))*100
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs*100
        
        t1CO_a[yy] = t1CO_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,7]))*100
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs*100
        t1CO_b[yy] = t1CO_b[yy] + 1*sum(L[[i]][[yy+1]][,6]>0)*100
        
        t2CO_a[yy] = t2CO_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,7]))*100
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs*100
        
        if(truen_num_blobs_ge4 > 0)
        {
          t1Hy_a[yy] = t1Hy_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 10] ))*100
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4*100
          t1Hy_b[yy] = t1Hy_b[yy] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11]) *100
          
          t2Hy_a[yy] = t2Hy_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 10] ))*100
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4*100
        }
      }
    }
    
    for(yy in 3:4)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        t2blobs_a[yy] = t2blobs_a[yy] + 0
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs*100
        
        t1blobs_a[yy] = t1blobs_a[yy] + 0
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs*100
        t1blobs_b[yy] = t1blobs_b[yy] + 0
        
        t2CO_a[yy] = t2CO_a[yy] + 0
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs*100
        
        t1CO_a[yy] = t1CO_a[yy] + 0
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs*100
        t1CO_b[yy] = t1CO_b[yy] + 0
        
        if(truen_num_blobs_ge4 > 0)
        {
          t2Hy_a[yy] = t2Hy_a[yy] + 0
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4*100
          t1Hy_b[yy] = t1Hy_b[yy] + 1*0
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4*100
        }
      } else {
        L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
        
        if(0%in%L[[i]][[yy+1]][,3])
        {
          L[[i]][[yy+1]] <- merge_zero_rows(L[[i]][[yy+1]])
        }
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,4]) * as.numeric(L[[i]][[yy+1]][,9]>0) )*100
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs*100
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0) *100
        
        t2blobs_a[yy] = t2blobs_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,4]) * as.numeric(L[[i]][[yy+1]][,9]>0) )*100
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs*100
        
        t1CO_a[yy] = t1CO_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,7]) )*100
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs*100
        t1CO_b[yy] = t1CO_b[yy] + 1*sum(L[[i]][[yy+1]][,6]>0)*100
        
        t2CO_a[yy] = t2CO_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,7]) )*100
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs*100
        
        if(truen_num_blobs_ge4 > 0)
        {
          t1Hy_a[yy] = t1Hy_a[yy] +  sum( as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 10]) )*100
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4*100
          t1Hy_b[yy] = t1Hy_b[yy] + 1*sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 9]>0)*100
          
          t2Hy_a[yy] = t2Hy_a[yy] +   sum( as.numeric(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 10]) )*100
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4*100
        }
      }
    }
  }
  
  idx = c(5,1,3,4,2)
  
  t1blobs  = (t1blobs_b[idx]+rep(sum(NN[,1]),5)-2*t1blobs_a[idx]) / (t1blobs_b[idx]+rep(sum(NN[,1]),5) -t1blobs_a[idx])
  t1CO     = (t1CO_b[idx]+   rep(sum(NN[,1]),5) -2*t1CO_a[idx]) /   (t1CO_b[idx]+rep(sum(NN[,1]),5) -t1CO_a[idx])
  t1Hy     = (t1Hy_b[idx]+   rep(sum(NN[,2]),5) -2*t1Hy_a[idx]) /   (t1Hy_b[idx]+rep(sum(NN[,2]),5) -t1Hy_a[idx])
  t1Sp     = (t1Sp_a/20)[idx]
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
  
  list(t1blobs = t1blobs, t1CO = t1CO, t1Hy = t1Hy, t1Sp = t1Sp, t1Mn = t1Mn)
}

 
compute_unrooted_metrics <- function(NN, LU, MU) {
  
  t1blobs_a = rep(0,3); t1blobs_b = rep(0,3); t1blobs_c = rep(0,3)
  t1CO_a    = rep(0,3); t1CO_b    = rep(0,3); t1CO_c    = rep(0,3)
  t1Hy_a    = rep(0,3); t1Hy_b    = rep(0,3); t1Hy_c    = rep(0,3)
  t1Sp      = rep(0,3)
  t2blobs_a = rep(0,3); t2blobs_b = rep(0,3); t2blobs_c = rep(0,3)
  t2CO_a    = rep(0,3); t2CO_b    = rep(0,3); t2CO_c    = rep(0,3)
  t2Hy_a    = rep(0,3); t2Hy_b    = rep(0,3); t2Hy_c    = rep(0,3)
  t2Sp      = rep(0,3)
  
  t1Sp_a = rep(0,3); t1Sp_b = rep(0,3); t1Sp_c = rep(0,3)
  t2Sp_a = rep(0,3); t2Sp_b = rep(0,3); t2Sp_c = rep(0,3)
  
  t1Sp_a[1]=sum(MU[,19]); t1Sp_b[1]=sum(1*MU[,1]); t1Sp_c = rep(sum(1*MU[,1]), 3)
  t1Sp_a[2]=sum(MU[,18]); t1Sp_b[2]=sum(1*MU[,2])
  t1Sp_a[3]=sum(MU[,17]); t1Sp_b[3]=sum(1*MU[,4])
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(LU[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(LU[[i]][[1]][LU[[i]][[1]][, 2] > 4, 8])
    
    yy=0
    if(LU[[i]][[yy+1]][1,1]=='NA'){
      flag=1
    } else {
      LU[[i]][[yy+1]][is.na(LU[[i]][[yy+1]])] <- 0
      
      t1blobs_a[yy+3] = t1blobs_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][,4]))*100
      t1blobs_c[yy+3] = t1blobs_c[yy+3] + 1*truen_num_blobs*100
      t1blobs_b[yy+3] = t1blobs_b[yy+3] + 1*sum(LU[[i]][[yy+1]][,3]>0)*100
      
      t2blobs_a[yy+3] = t2blobs_a[yy+3] +  sum(as.numeric(LU[[i]][[yy+1]][,4]))*100
      t2blobs_c[yy+3] = t2blobs_c[yy+3] + 1*truen_num_blobs*100
      
      t1CO_a[yy+3] = t1CO_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][,6]))*100
      t1CO_c[yy+3] = t1CO_c[yy+3] + 1*truen_num_blobs*100
      t1CO_b[yy+3] = t1CO_b[yy+3] + 1*sum(LU[[i]][[yy+1]][,6]>0)*100
      
      t2CO_a[yy+3] = t2CO_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][,6]))*100
      t2CO_c[yy+3] = t2CO_c[yy+3] + 1*truen_num_blobs*100
      
      if(truen_num_blobs_ge4 > 0)
      {
        t1Hy_a[yy+3] = t1Hy_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 8] ))*100
        t1Hy_c[yy+3] = t1Hy_c[yy+3] + 1*truen_num_blobs_ge4*100
        t1Hy_b[yy+3] = t1Hy_b[yy+3] + 1*length(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 7]) *100
        
        t2Hy_a[yy+3] = t2Hy_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 8] ))*100
        t2Hy_c[yy+3] = t2Hy_c[yy+3] + 1*truen_num_blobs_ge4*100
      }
    }
    
    for(yy in 1:2)
    {
      if(LU[[i]][[yy+1]][1,1]=='NA'){
        flag=1
      } else {
        LU[[i]][[yy+1]][is.na(LU[[i]][[yy+1]])] <- 0
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][,4]))*100
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs*100
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(LU[[i]][[yy+1]][,9]>0)*100
        
        t2blobs_a[yy] = t2blobs_a[yy] +  sum(as.numeric(LU[[i]][[yy+1]][,4]))*100
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs*100
        
        t1CO_a[yy] = t1CO_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][,7]))*100
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs*100
        t1CO_b[yy] = t1CO_b[yy] + 1*sum(LU[[i]][[yy+1]][,6]>0)*100
        
        t2CO_a[yy] = t2CO_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][,7]))*100
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs*100
        
        if(truen_num_blobs_ge4 > 0)
        {
          t1Hy_a[yy] = t1Hy_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 10] ))*100
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4*100
          t1Hy_b[yy] = t1Hy_b[yy] + 1*length(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 11]) *100
          
          t2Hy_a[yy] = t2Hy_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 10] ))*100
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4*100
        }
      }
    }
  }
  
  idx = c(3,1,2)
  
  
  
  t1blobs  = (t1blobs_b[idx]+rep(sum(NN[,1]),3) -2*t1blobs_a[idx]) / (t1blobs_b[idx]+rep(sum(NN[,1]),3) -t1blobs_a[idx])
  t1CO     = (t1CO_b[idx]   +rep(sum(NN[,1]),3) -2*t1CO_a[idx]) / (t1CO_b[idx]+rep(sum(NN[,1]),3) -t1CO_a[idx])
  t1Hy     = (t1Hy_b[idx]   +rep(sum(NN[,2]),3) -2*t1Hy_a[idx]) / (t1Hy_b[idx]+rep(sum(NN[,2]),3) -t1Hy_a[idx])
  t1Sp    = (t1Sp_a/20)[idx]
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
  
  list(t1blobs = t1blobs, t1CO = t1CO, t1Hy = t1Hy, t1Sp = t1Sp, t1Mn = t1Mn)
}

 
metrics_names <- c("t1blobs","t1CO","t1Hy","t1Sp","t1Mn")

rep_rooted <- vector("list", 2)
rep_unroot <- vector("list", 2)
for (kk in 1:2) {
  rep_rooted[[kk]] <- setNames(lapply(metrics_names, function(nm) matrix(NA_real_, nrow = n_reps, ncol = 5)), metrics_names)
  rep_unroot[[kk]] <- setNames(lapply(metrics_names, function(nm) matrix(NA_real_, nrow = n_reps, ncol = 3)), metrics_names)
}

for (rep_i in 1:n_reps) {
  
  message("=== Replicate ", rep_i, " ===")
  
  dir_rootedfull    <- file.path(sim_out, paste0("rootedfull", rep_i))
  dir_rootedthird   <- file.path(sim_out, paste0("rootedthird", rep_i))
  dir_unrootedfull  <- file.path(sim_out, paste0("unrootedfull", rep_i))
  dir_unrootedthird <- file.path(sim_out, paste0("unrootedthird", rep_i))
  
  message("  computing bootstrap blob counts (NB matrices)...")
  NBfull   <- compute_NB_matrix(dir_rootedfull,    function(xp) paste0("rooted_network_newick_concatalignRealisticSim_400gts_Net", xp, ".txt"))
  NBthird  <- compute_NB_matrix(dir_rootedthird,   function(xp) paste0("rooted_network_newick_concatalignRealisticSim_400gts_Net", xp, "_third.txt"))
  NBUfull  <- compute_NB_matrix(dir_unrootedfull,  function(xp) paste0("unrooted_network_newick_concatalignRealisticSim_400gts_Net", xp, ".txt"))
  NBUthird <- compute_NB_matrix(dir_unrootedthird, function(xp) paste0("unrooted_network_newick_concatalignRealisticSim_400gts_Net", xp, "_third.txt"))
  
  message("  loading Blob/Edge analysis Rdata...")
  load(file.path(dir_unrootedfull,  "BlobAnalysis_full_nooutgroup.Rdata"));  result1 <- result
  load(file.path(dir_unrootedthird, "BlobAnalysis_third_nooutgroup.Rdata")); result2 <- result
  load(file.path(dir_rootedfull,    "BlobAnalysis_full_outgroup.Rdata"));    result3 <- result
  load(file.path(dir_rootedthird,   "BlobAnalysis_third_outgroup.Rdata"));   result4 <- result
  
  load(file.path(dir_rootedfull,    "EdgeAnalysis_full_outgroup_jaccard.Rdata"))    # -> MM1
  load(file.path(dir_rootedthird,   "EdgeAnalysis_third_outgroup_jaccard.Rdata"))   # -> MM2
  load(file.path(dir_unrootedfull,  "EdgeAnalysis_full_nooutgroup_jaccard.Rdata"))  # -> MM3
  load(file.path(dir_unrootedthird, "EdgeAnalysis_third_nooutgroup_jaccard.Rdata")) # -> MM4
  
  for (kk in 1:2) {
    
    NN <- list(NBfull, NBthird)[[kk]]
    L  <- list(result3, result4)[[kk]]
    M  <- list(MM1, MM2)[[kk]]
    
    res_rooted <- compute_rooted_metrics(NN, L, M)
    for (nm in metrics_names) rep_rooted[[kk]][[nm]][rep_i, ] <- res_rooted[[nm]]
    
    NNu <- list(NBUfull, NBUthird)[[kk]]
    LU  <- list(result1, result2)[[kk]]
    MU  <- list(MM3, MM4)[[kk]]
    
    res_unroot <- compute_unrooted_metrics(NNu, LU, MU)
    for (nm in metrics_names) rep_unroot[[kk]][[nm]][rep_i, ] <- res_unroot[[nm]]
  }
}

 
avg_rooted <- vector("list", 2)
avg_unroot <- vector("list", 2)
for (kk in 1:2) {
  avg_rooted[[kk]] <- lapply(rep_rooted[[kk]], function(m) colMeans(m, na.rm = TRUE))
  avg_unroot[[kk]] <- lapply(rep_unroot[[kk]], function(m) colMeans(m, na.rm = TRUE))
}

save(rep_rooted, rep_unroot, avg_rooted, avg_unroot,
     file = file.path(sim_out, "AveragedInputDissimilarity_10reps.Rdata"))

 
plotsroot  <- list()
plotsunroot <- list()

for (kk in 1:2) {
  
  t1blobs <- avg_rooted[[kk]]$t1blobs
  t1CO    <- avg_rooted[[kk]]$t1CO
  t1Hy    <- avg_rooted[[kk]]$t1Hy
  t1Sp    <- avg_rooted[[kk]]$t1Sp
  t1Mn    <- avg_rooted[[kk]]$t1Mn
  
  df1 <- data.frame(
    value = c(t1blobs, t1CO, t1Hy, t1Sp, t1Mn),
    category = rep(category_rooted, 5),
    measure  = rep(measure_levels, each = 5)
  )
  df1$category <- factor(df1$category, levels = category_rooted)
  df1$measure  <- factor(df1$measure,  levels = measure_levels)
  
  p1 <- ggplot(df1, aes(x = measure, y = value, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols_rooted) +
    labs(
      title = paste0(c("(a) ","(c) ")[kk], "Jaccard Dissimilarity to Input networks - Rooted - ",
                     c("Full","Third")[kk], " (avg. of ", n_reps, " reps)"),
      x = "", y = "", fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", plot.title = element_text(hjust = 0.5)
    ) +
    big_text +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1))
  
  plotsroot[[kk]] <- p1
}

for (kk in 1:2) {
  
  t1blobs <- avg_unroot[[kk]]$t1blobs
  t1CO    <- avg_unroot[[kk]]$t1CO
  t1Hy    <- avg_unroot[[kk]]$t1Hy
  t1Sp    <- avg_unroot[[kk]]$t1Sp
  t1Mn    <- avg_unroot[[kk]]$t1Mn
  
  df1 <- data.frame(
    value = c(t1blobs, t1CO, t1Hy, t1Sp, t1Mn),
    category = rep(category_unroot, 5),
    measure  = rep(measure_levels, each = 3)
  )
  df1$category <- factor(df1$category, levels = category_unroot)
  df1$measure  <- factor(df1$measure,  levels = measure_levels)
  
  p1 <- ggplot(df1, aes(x = measure, y = value, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols_unroot) +
    labs(
      title = paste0(c("(b) ","(d) ")[kk], "Jaccard Dissimilarity to Input networks - Unrooted - ",
                     c("Full","Third")[kk], " (avg. of ", n_reps, " reps)"),
      x = "", y = "", fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", plot.title = element_text(hjust = 0.5)
    ) +
    big_text +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1))
  
  plotsunroot[[kk]] <- p1
}

combined_plot <- grid.arrange(
  grobs = c(plotsroot[[1]], plotsunroot[[1]]),
  ncol = 2
)
ggsave(
  file.path(sim_out, "DistFull-Average_jaccard.pdf"),
  plot = combined_plot,
  width = 19,
  height = 6.5,
  units = "in"
)

combined_plot <- grid.arrange(
  grobs = c(plotsroot[[2]], plotsunroot[[2]]),
  ncol = 2
)
ggsave(
  file.path(sim_out, "DistThird-Average_jaccard.pdf"),
  plot = combined_plot,
  width = 19,
  height = 6.5,
  units = "in"
)

message("Done. Averaged plots and AveragedInputDissimilarity_10reps.Rdata written to: ", sim_out)