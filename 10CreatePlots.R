###########################################################################################
#################### Create plots for manuscript (PRec, Recall, F1) ######################
###########################################################################################

### This script was written by a human and optimized using Claude LLM



library(openxlsx)
library(ggplot2)
library(patchwork)

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"
sim_out  <- file.path(base_dir, "sim_out")

n_reps <- 10

category_levels <- c("Squirrel","HMS","HWC","L1N-majority","L1N-unrooted","L1N")
measure_levels  <- c("Blobs","Cycles","Hybrids B>4","Splits","Mean")

cols <- c(
  "gray70",
  "#A8DADC",
  "#274C77",
  "#A8D5A2",
  "#4CAF50",
  "#1B5E20"
)

  
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

 
compute_kk_metrics <- function(L, LU, M, MU) {
  
  t1blobs_a = rep(0,6); t1blobs_b = rep(0,6); t1blobs_c = rep(0,6)
  t1CO_a    = rep(0,6); t1CO_b    = rep(0,6); t1CO_c    = rep(0,6)
  t1Hy_a    = rep(0,6); t1Hy_b    = rep(0,6); t1Hy_c    = rep(0,6)
  t1Sp      = rep(0,6)
  t2blobs_a = rep(0,6); t2blobs_b = rep(0,6); t2blobs_c = rep(0,6)
  t2CO_a    = rep(0,6); t2CO_b    = rep(0,6); t2CO_c    = rep(0,6)
  t2Hy_a    = rep(0,6); t2Hy_b    = rep(0,6); t2Hy_c    = rep(0,6)
  t2Sp      = rep(0,6)
  
  t1Sp_a = rep(0,6); t1Sp_b = rep(0,6); t1Sp_c = rep(0,6)
  t2Sp_a = rep(0,6); t2Sp_b = rep(0,6); t2Sp_c = rep(0,6)
  
  bg4  = 0
  sipo = c()
  
  t1Sp_a[1]=sum(M[,15]); t1Sp_c[1]=sum(1*M[,1]); t1Sp_b[1]=sum(1*M[,3])
  t2Sp_a[1]=sum(M[,15]); t2Sp_c[1]=sum(1*M[,1]); t2Sp_b[1]=sum(1*M[,3])
  
  t1Sp_a[2]=sum(M[,12]); t1Sp_c[2]=sum(1*M[,1]); t1Sp_b[2]=sum(1*M[,2])
  t2Sp_a[2]=sum(M[,12]); t2Sp_c[2]=sum(1*M[,1]); t2Sp_b[2]=sum(1*M[,2])
  
  t1Sp_a[3]=sum(M[,18]); t1Sp_c[3]=sum(1*M[,1]); t1Sp_b[3]=sum(1*M[,4])
  t2Sp_a[3]=sum(M[,18]); t2Sp_c[3]=sum(1*M[,1]); t2Sp_b[3]=sum(1*M[,4])
  
  t1Sp_a[5]=sum(M[,24]); t1Sp_c[4]=sum(1*M[,1]); t1Sp_b[5]=sum(1*M[,6])
  t2Sp_a[5]=sum(M[,24]); t2Sp_c[4]=sum(1*M[,1]); t2Sp_b[5]=sum(1*M[,6])
  
  t1Sp_a[4]=sum(M[,21]); t1Sp_c[5]=sum(1*M[,1]); t1Sp_b[4]=sum(1*M[,5])
  t2Sp_a[4]=sum(M[,21]); t2Sp_c[5]=sum(1*M[,1]); t2Sp_b[4]=sum(1*M[,5])
  
  t1Sp_a[6]=sum(MU[,8]); t1Sp_c[6]=sum(1*MU[,1]); t1Sp_b[6]=sum(1*MU[,2])
  t2Sp_a[6]=sum(MU[,8]); t2Sp_c[6]=sum(1*MU[,1]); t2Sp_b[6]=sum(1*MU[,2])
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(L[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(L[[i]][[1]][L[[i]][[1]][, 2] > 4, 8])
    
    for(yy in 1:2)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        t2blobs_a[yy] = t2blobs_a[yy] + 0
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs
        if(yy==3){sipo =c(sipo, 0)}
        t1blobs_a[yy] = t1blobs_a[yy] + 0
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs
        t1blobs_b[yy] = t1blobs_b[yy] + 0
        t2CO_a[yy] = t2CO_a[yy] + 0
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs
        t1CO_a[yy] = t1CO_a[yy] + 0
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs
        t1CO_b[yy] = t1CO_b[yy] + 0
        if(truen_num_blobs_ge4 > 0)
        {
          t2Hy_a[yy] = t2Hy_a[yy] + 0
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4
          t1Hy_b[yy] = t1Hy_b[yy] + 1*0
          bg4 = bg4+1
        }
      } else {
        L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(L[[i]][[yy+1]][,5]==1)
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0)
        
        t2blobs_a[yy] = t2blobs_a[yy] +   sum(L[[i]][[yy+1]][,5]==1)
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs
        
        t1CO_a[yy] = t1CO_a[yy] +   sum(L[[i]][[yy+1]][,8]==1)
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs
        t1CO_b[yy] = t1CO_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0)
        
        t2CO_a[yy] = t2CO_a[yy] +   sum(L[[i]][[yy+1]][,8]==1)
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs
        
        if(truen_num_blobs_ge4 > 0)
        {
          t1Hy_a[yy] = t1Hy_a[yy] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11]==1)
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4
          t1Hy_b[yy] = t1Hy_b[yy] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11])
          
          t2Hy_a[yy] = t2Hy_a[yy] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11]==1)
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4
          bg4 = bg4+1
        }
      }
    }
    
    for(yy in 3:4)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        t2blobs_a[yy] = t2blobs_a[yy] + 0
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs
        if(yy==3){sipo =c(sipo, 0)}
        t1blobs_a[yy] = t1blobs_a[yy] + 0
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs
        t1blobs_b[yy] = t1blobs_b[yy] + 0
        t2CO_a[yy] = t2CO_a[yy] + 0
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs
        t1CO_a[yy] = t1CO_a[yy] + 0
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs
        t1CO_b[yy] = t1CO_b[yy] + 0
        if(truen_num_blobs_ge4 > 0)
        {
          t2Hy_a[yy] = t2Hy_a[yy] + 0
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4
          t1Hy_b[yy] = t1Hy_b[yy] + 1*0
          bg4 = bg4+1
        }
      } else {
        L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
        
        if(0%in%L[[i]][[yy+1]][,3])
        {
          
          L[[i]][[yy+1]] <- merge_zero_rows(L[[i]][[yy+1]])
        }
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(L[[i]][[yy+1]][,5]==1 & L[[i]][[yy+1]][,9]>0)
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0)
        if(yy==3){sipo =c(sipo, sum(L[[i]][[yy+1]][,5]==1 & L[[i]][[yy+1]][,9]>0))}
        t2blobs_a[yy] = t2blobs_a[yy] +   sum(L[[i]][[yy+1]][,5]==1 & L[[i]][[yy+1]][,9]>0)
        t2blobs_c[yy] = t2blobs_c[yy] + 1*truen_num_blobs
        
        t1CO_a[yy] = t1CO_a[yy] +   sum(L[[i]][[yy+1]][,8]==1)
        t1CO_c[yy] = t1CO_c[yy] + 1*truen_num_blobs
        t1CO_b[yy] = t1CO_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0)
        
        t2CO_a[yy] = t2CO_a[yy] +   sum(L[[i]][[yy+1]][,8]==1)
        t2CO_c[yy] = t2CO_c[yy] + 1*truen_num_blobs
        
        if(truen_num_blobs_ge4 > 0)
        {
          t1Hy_a[yy] = t1Hy_a[yy] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11]==1)
          t1Hy_c[yy] = t1Hy_c[yy] + 1*truen_num_blobs_ge4
          t1Hy_b[yy] = t1Hy_b[yy] + 1*sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 9]>0)
          
          t2Hy_a[yy] = t2Hy_a[yy] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4, 11]==1)
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4
          bg4 = bg4+1
        }
      }
    }
    
    yy = 2
    
    if(L[[i]][[yy+1]][1,1]=='NA'){
      t2blobs_a[5] = t2blobs_a[5] + 0
      t2blobs_c[5] = t2blobs_c[5] + 1*truen_num_blobs
      if(yy==3){sipo =c(sipo, 0)}
      t1blobs_a[5] = t1blobs_a[5] + 0
      t1blobs_c[5] = t1blobs_c[5] + 1*truen_num_blobs
      t1blobs_b[5] = t1blobs_b[5] + 0
      t2CO_a[5] = t2CO_a[5] + 0
      t2CO_c[5] = t2CO_c[5] + 1*truen_num_blobs
      t1CO_a[5] = t1CO_a[5] + 0
      t1CO_c[5] = t1CO_c[5] + 1*truen_num_blobs
      t1CO_b[5] = t1CO_b[5] + 0
      if(truen_num_blobs_ge4 > 0)
      {
        t2Hy_a[5] = t2Hy_a[5] + 0
        t2Hy_c[5] = t2Hy_c[5] + 1*truen_num_blobs_ge4
        t1Hy_b[5] = t1Hy_b[5] + 1*0
        bg4 = bg4+1
      }
    } else {
      L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
      
      t1blobs_a[5] = t1blobs_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,4]>0.5, 5]==1)
      t1blobs_c[5] = t1blobs_c[5] + 1*truen_num_blobs
      t1blobs_b[5] = t1blobs_b[5] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,4]>0.5, 5])
      
      t2blobs_a[5] = t2blobs_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,4]>0.5, 5]==1)
      t2blobs_c[5] = t2blobs_c[5] + 1*truen_num_blobs
      
      t1CO_a[5] = t1CO_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,7]>0.5, 8]==1)
      t1CO_c[5] = t1CO_c[5] + 1*truen_num_blobs
      t1CO_b[5] = t1CO_b[5] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,7]>0.5, 5])
      
      t2CO_a[5] = t2CO_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,7]>0.5, 8]==1)
      t2CO_c[5] = t2CO_c[5] + 1*truen_num_blobs
      
      if(truen_num_blobs_ge4 > 0)
      {
        t1Hy_a[5] = t1Hy_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4 & L[[i]][[yy+1]][,10]>0.5, 11]==1)
        t1Hy_c[5] = t1Hy_c[5] + 1*truen_num_blobs_ge4
        t1Hy_b[5] = t1Hy_b[5] + 1*length(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4 & L[[i]][[yy+1]][,10]>0.5, 11])
        
        t2Hy_a[5] = t2Hy_a[5] +   sum(L[[i]][[yy+1]][L[[i]][[yy+1]][,2]>4 & L[[i]][[yy+1]][,10]>0.5, 11]==1)
        t2Hy_c[5] = t2Hy_c[5] + 1*truen_num_blobs_ge4
      }
    }
    
    yy = 2
    
    if(LU[[i]][[yy+1]][1,1]=='NA'){
      t2blobs_a[6] = t2blobs_a[6] + 0
      t2blobs_c[6] = t2blobs_c[6] + 1*truen_num_blobs
      if(yy==3){sipo =c(sipo, 0)}
      t1blobs_a[6] = t1blobs_a[6] + 0
      t1blobs_c[6] = t1blobs_c[6] + 1*truen_num_blobs
      t1blobs_b[6] = t1blobs_b[6] + 0
      t2CO_a[6] = t2CO_a[6] + 0
      t2CO_c[6] = t2CO_c[6] + 1*truen_num_blobs
      t1CO_a[6] = t1CO_a[6] + 0
      t1CO_c[6] = t1CO_c[6] + 1*truen_num_blobs
      t1CO_b[6] = t1CO_b[6] + 0
      if(truen_num_blobs_ge4 > 0)
      {
        t2Hy_a[6] = t2Hy_a[6] + 0
        t2Hy_c[6] = t2Hy_c[6] + 1*truen_num_blobs_ge4
        t1Hy_b[6] = t1Hy_b[6] + 1*0
        bg4 = bg4+1
      }
    } else {
      LU[[i]][[yy+1]][is.na(LU[[i]][[yy+1]])] <- 0
      
      t1blobs_a[6] = t1blobs_a[6] +   sum(LU[[i]][[yy+1]][,5]==1)
      t1blobs_c[6] = t1blobs_c[6] + 1*truen_num_blobs
      t1blobs_b[6] = t1blobs_b[6] + 1*length(LU[[i]][[yy+1]][,1])
      
      t2blobs_a[6] = t2blobs_a[6] +   sum(LU[[i]][[yy+1]][,5]==1)
      t2blobs_c[6] = t2blobs_c[6] + 1*truen_num_blobs
      
      t1CO_a[6] = t1CO_a[6] +   sum(LU[[i]][[yy+1]][,8]==1)
      t1CO_c[6] = t1CO_c[6] + 1*truen_num_blobs
      t1CO_b[6] = t1CO_b[6] + 1*length(LU[[i]][[yy+1]][,1])
      
      t2CO_a[6] = t2CO_a[6] +   sum(LU[[i]][[yy+1]][,8]==1)
      t2CO_c[6] = t2CO_c[6] + 1*truen_num_blobs
      
      if(truen_num_blobs_ge4 > 0)
      {
        t1Hy_a[6] = t1Hy_a[6] +   sum(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 11]==1)
        t1Hy_c[6] = t1Hy_c[6] + 1*truen_num_blobs_ge4
        t1Hy_b[6] = t1Hy_b[6] + 1*length(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 9])
        
        t2Hy_a[6] = t2Hy_a[6] +   sum(LU[[i]][[yy+1]][LU[[i]][[yy+1]][,2]>4, 11]==1)
        t2Hy_c[6] = t2Hy_c[6] + 1*truen_num_blobs_ge4
      }
    }
  }
  
  idx = c(1,3,4,5,6,2)
  
  t1blobs = (t1blobs_a / t1blobs_b)[idx]
  t2blobs = (t2blobs_a / t2blobs_c)[idx]
  
  t1CO = (t1CO_a / t1CO_b)[idx]
  t2CO = (t2CO_a / t2CO_c)[idx]
  
  t1Hy = (t1Hy_a / t1Hy_b)[idx]
  t2Hy = (t2Hy_a / t2Hy_c)[idx]
  
  t1Sp = (t1Sp_a / t1Sp_b)[idx]
  t2Sp = (t2Sp_a / t2Sp_c)[idx]
  
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
  t2Mn = (t2blobs + t2CO + t2Hy + t2Sp) / 4
  
  dblobs = (2 * (t1blobs) * (t2blobs)) / ((t1blobs) + (t2blobs))
  dCO    =  2 * (t1CO) * (t2CO) / ((t1CO) + (t2CO))
  dHy    =  2 * (t1Hy) * (t2Hy) / ((t1Hy) + (t2Hy))
  dSp    =   2 * (t1Sp) * (t2Sp) / ((t1Sp) + (t2Sp))
  dMn    = (dblobs + dCO + dHy + dSp) / 4
  
  drfblobs = (t1blobs_b[idx]+t1blobs_c[idx] -2*t1blobs_a[idx]) / (t1blobs_b[idx]+t1blobs_c[idx] -t1blobs_a[idx])
  drfCO    = (t1CO_b[idx]+t1CO_c[idx] -2*t1CO_a[idx]) / (t1CO_b[idx]+t1CO_c[idx] -t1CO_a[idx])
  drfHy    = (t1Hy_b[idx]+t1Hy_c[idx] -2*t1Hy_a[idx]) / (t1Hy_b[idx]+t1Hy_c[idx] -t1Hy_a[idx])
  drfSp    = (t1Sp_b[idx]+t1Sp_c[idx] -2*t1Sp_a[idx]) / (t1Sp_b[idx]+t1Sp_c[idx] -t1Sp_a[idx])
  drfMn    = (drfblobs + drfCO + drfHy + drfSp) / 4
  
  list(
    t1blobs = t1blobs, t2blobs = t2blobs,
    t1CO = t1CO,       t2CO = t2CO,
    t1Hy = t1Hy,       t2Hy = t2Hy,
    t1Sp = t1Sp,       t2Sp = t2Sp,
    t1Mn = t1Mn,       t2Mn = t2Mn,
    dblobs = dblobs,   dCO = dCO,   dHy = dHy,   dSp = dSp,   dMn = dMn,
    drfblobs = drfblobs, drfCO = drfCO, drfHy = drfHy, drfSp = drfSp, drfMn = drfMn
  )
}

 
metrics_names <- c("t1blobs","t2blobs","t1CO","t2CO","t1Hy","t2Hy",
                   "t1Sp","t2Sp","t1Mn","t2Mn",
                   "dblobs","dCO","dHy","dSp","dMn",
                   "drfblobs","drfCO","drfHy","drfSp","drfMn")

rep_metrics <- vector("list", 2)
for (kk in 1:2) {
  rep_metrics[[kk]] <- setNames(
    lapply(metrics_names, function(nm) matrix(NA_real_, nrow = n_reps, ncol = 6)),
    metrics_names
  )
}

for (rep_i in 1:n_reps) {
  
  message("=== Loading + computing replicate ", rep_i, " ===")
  
   
  load(file.path(sim_out, paste0("unrootedfull",  rep_i), "BlobAnalysis_full_nooutgroup.Rdata"));  result1 <- result
  load(file.path(sim_out, paste0("unrootedthird", rep_i), "BlobAnalysis_third_nooutgroup.Rdata")); result2 <- result
  load(file.path(sim_out, paste0("rootedfull",    rep_i), "BlobAnalysis_full_outgroup.Rdata"));    result3 <- result
  load(file.path(sim_out, paste0("rootedthird",   rep_i), "BlobAnalysis_third_outgroup.Rdata"));   result4 <- result
  
   
  load(file.path(sim_out, paste0("rootedfull",    rep_i), "EdgeAnalysis_full_outgroup.Rdata"))    # -> MM1
  load(file.path(sim_out, paste0("rootedthird",   rep_i), "EdgeAnalysis_third_outgroup.Rdata"))   # -> MM2
  load(file.path(sim_out, paste0("unrootedfull",  rep_i), "EdgeAnalysis_full_nooutgroup.Rdata"))  # -> MM3
  load(file.path(sim_out, paste0("unrootedthird", rep_i), "EdgeAnalysis_third_nooutgroup.Rdata")) # -> MM4
  
  for (kk in 1:2) {
    L  <- list(result3, result4)[[kk]]
    LU <- list(result1, result2)[[kk]]
    M  <- list(MM1, MM2)[[kk]]
    MU <- list(MM3, MM4)[[kk]]
    
    res <- compute_kk_metrics(L, LU, M, MU)
    
    for (nm in metrics_names) {
      rep_metrics[[kk]][[nm]][rep_i, ] <- res[[nm]]
    }
  }
}

 
avg_metrics <- vector("list", 2)
for (kk in 1:2) {
  avg_metrics[[kk]] <- lapply(rep_metrics[[kk]], function(m) colMeans(m, na.rm = TRUE))
}

save(rep_metrics, avg_metrics, file = file.path(sim_out, "AveragedMetrics_PR_10reps.Rdata"))

 
distanceplot <- list()

for (kk in 1:2) {
  
  t1blobs <- avg_metrics[[kk]]$t1blobs; t2blobs <- avg_metrics[[kk]]$t2blobs
  t1CO    <- avg_metrics[[kk]]$t1CO;    t2CO    <- avg_metrics[[kk]]$t2CO
  t1Hy    <- avg_metrics[[kk]]$t1Hy;    t2Hy    <- avg_metrics[[kk]]$t2Hy
  t1Sp    <- avg_metrics[[kk]]$t1Sp;    t2Sp    <- avg_metrics[[kk]]$t2Sp
  t1Mn    <- avg_metrics[[kk]]$t1Mn;    t2Mn    <- avg_metrics[[kk]]$t2Mn
  dblobs  <- avg_metrics[[kk]]$dblobs;  dCO     <- avg_metrics[[kk]]$dCO
  dHy     <- avg_metrics[[kk]]$dHy;     dSp     <- avg_metrics[[kk]]$dSp
  dMn     <- avg_metrics[[kk]]$dMn
  drfblobs <- avg_metrics[[kk]]$drfblobs; drfCO <- avg_metrics[[kk]]$drfCO
  drfHy    <- avg_metrics[[kk]]$drfHy;    drfSp <- avg_metrics[[kk]]$drfSp
  drfMn    <- avg_metrics[[kk]]$drfMn
  
  distanceplot[[kk]] <- c(drfblobs, drfCO, drfHy, drfSp, drfMn)
  
  df1 <- data.frame(
    value = c(t1blobs, t1CO, t1Hy, t1Sp, t1Mn),
    category = rep(category_levels, 5),
    measure  = rep(measure_levels, each = 6)
  )
  df1$category <- factor(df1$category, levels = category_levels)
  df1$measure  <- factor(df1$measure,  levels = measure_levels)
  
  df2 <- data.frame(
    value = c(t2blobs, t2CO, t2Hy, t2Sp, t2Mn),
    category = rep(category_levels, 5),
    measure  = rep(measure_levels, each = 6)
  )
  df2$category <- factor(df2$category, levels = category_levels)
  df2$measure  <- factor(df2$measure,  levels = measure_levels)
  
  dfd <- data.frame(
    value = c(dblobs, dCO, dHy, dSp, dMn),
    category = rep(category_levels, 5),
    measure  = rep(measure_levels, each = 6)
  )
  dfd$category <- factor(dfd$category, levels = category_levels)
  dfd$measure  <- factor(dfd$measure,  levels = measure_levels)
  
  p1 <- ggplot(df1, aes(x = measure, y = value, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0("Precision - ", c("Full length","Third of each gene")[kk],
                     " (avg. of ", n_reps, " reps)"),
      x = "", y = "", fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", plot.title = element_text(hjust = 0.5)
    ) +
    big_text +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = "none")
  
  p2 <- ggplot(df2, aes(x = measure, y = value, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0("Recall - ", c("Full length","Third of each gene")[kk],
                     " (avg. of ", n_reps, " reps)"),
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
  
  p3 <- ggplot(dfd, aes(x = measure, y = value, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = bquote(F[1]~score~"-"~.(c("Full length","Third of each gene")[kk])~
                       "(avg. of"~.(n_reps)~"reps)"),
      x = "", y = "", fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", plot.title = element_text(hjust = 0.5)
    ) +
    big_text +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = "none")
  
  print(p1 + p2 + p3)
  
  ggsave(
    file.path(sim_out, paste0(c("Full","Third")[kk], "-Average-data.pdf")),
    plot = p1 + p2 + p3,
    width = 19,
    height = 6.5,
    units = "in"
  )
}

 
dfd1 <- data.frame(
  value = distanceplot[[1]],
  category = rep(category_levels, 5),
  measure  = rep(measure_levels, each = 6)
)
dfd1$category <- factor(dfd1$category, levels = category_levels)
dfd1$measure  <- factor(dfd1$measure,  levels = measure_levels)

dfd2 <- data.frame(
  value = distanceplot[[2]],
  category = rep(category_levels, 5),
  measure  = rep(measure_levels, each = 6)
)
dfd2$category <- factor(dfd2$category, levels = category_levels)
dfd2$measure  <- factor(dfd2$measure,  levels = measure_levels)

pp1 <- ggplot(dfd1, aes(x = measure, y = value, fill = category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cols) +
  labs(
    title = paste0("(a) Jaccard Dissimilarity - ", c("Full length","Third of each gene")[1],
                   " (avg. of ", n_reps, " reps)"),
    x = "", y = "", fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.box = "horizontal", plot.title = element_text(hjust = 0.5)
  ) +
  big_text +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(fill = "none")

pp2 <- ggplot(dfd2, aes(x = measure, y = value, fill = category)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cols) +
  labs(
    title = paste0("(b) Jaccard Dissimilarity - ", c("Full length","Third of each gene")[2],
                   " (avg. of ", n_reps, " reps)"),
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

cp1 <- (pp1 + pp2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  ) &
  big_text

ggsave(
  file.path(sim_out, "DisttoTrue--Average.pdf"),
  plot = cp1,
  width = 19,
  height = 6.5,
  units = "in"
)
message("Done. Averaged plots and AveragedMetrics_PR_10reps.Rdata written to: ", sim_out)