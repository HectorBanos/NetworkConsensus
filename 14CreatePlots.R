library(openxlsx)
library(ggplot2)
library(patchwork)

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2"
location <- file.path(base_dir, "sim_out")
setwd(location)
file1=paste0(location,'/BlobAnalysis_full_nooutgroup.Rdata')
file2=paste0(location,'/BlobAnalysis_third_nooutgroup.Rdata')
file3=paste0(location,'/BlobAnalysis_full_outgroup.Rdata')
file4=paste0(location,'/BlobAnalysis_third_outgroup.Rdata')
load(file1)
result1=result
load(file2)
result2=result
load(file3)
result3=result
load(file4)
result4=result
distanceplot=list()
merge_zero_rows <- function(df) {
  score_cols <- setdiff(names(df), c("canon", "n_partitions"))
  first_col <- score_cols[1]
  
  repeat {
    zero_rows <- which(df[[first_col]] == 0)
    if (length(zero_rows) == 0) break
    
    i <- zero_rows[which.min(df$n_partitions[zero_rows])]
    s <- df$n_partitions[i]
    
    # Columns where row i has NON-zero values (these will be donated)
    nonzero_cols_i <- score_cols[as.numeric(df[i, score_cols]) != 0]
    
    if (length(nonzero_cols_i) == 0) break
    
    # Target: n_partitions > s, AND all of row i's non-zero cols are zero in target
    candidates <- which(
      df$n_partitions > s &
        apply(df[, nonzero_cols_i, drop = FALSE], 1,
              function(r) all(as.numeric(r) == 0))
    )
    
    if (length(candidates) == 0) break
    
    j <- candidates[which.min(df$n_partitions[candidates])]
    
    # Fill j's zeros with i's non-zero values
    for (col in nonzero_cols_i) {
      df[j, col] <- df[i, col]
    }
    
    df <- df[-i, ]
  }
  
  rownames(df) <- NULL
  return(df)
}

file3=paste0(location,'/EdgeAnalysis_full_nooutgroup.Rdata')
file4=paste0(location,'/EdgeAnalysis_third_nooutgroup.Rdata')
file1=paste0(location,'/EdgeAnalysis_full_outgroup.Rdata')
file2=paste0(location,'/EdgeAnalysis_third_outgroup.Rdata')

load(file3)#M3
load(file4)#M4
load(file1)#M1
load(file2)#M2

 sipo=list()
kk=1
for(kk in 1:2)
{
  L =list(result3,result4)[[kk]]
  LU=list(result1,result2)[[kk]]
  
  M=list(MM1,MM2)[[kk]]
  MU=list(MM3,MM4)[[kk]]
  
 
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
  
  bg4 = 0
  

  # t1Sp[1]=sum(M[,17])/sum(1*M[,4-1]); t2Sp[1]=sum(M[,16])/sum(1*M[,1])
  # t1Sp[2]=sum(M[,14])/sum(1*M[,3-1]); t2Sp[2]=sum(M[,13])/sum(1*M[,1])  
  # t1Sp[3]=sum(M[,20])/sum(1*M[,5-1]); t2Sp[3]=sum(M[,19])/sum(1*M[,1])
  # t1Sp[4]=sum(M[,26])/sum(1*M[,7-1]); t2Sp[4]=sum(M[,25])/sum(1*M[,1])
  # t1Sp[6]=sum(M[,23])/sum(1*M[,6-1]); t2Sp[6]=sum(M[,22])/sum(1*M[,1])
  # t1Sp[5]=sum(MU[,10])/sum(1*MU[,3-1]); t2Sp[5]=sum(MU[,9])/sum(1*MU[,1])
  
  t1Sp_a[1]=sum(M[,15]);
  t1Sp_c[1]=sum(1*M[,1]);
  t1Sp_b[1]=sum(1*M[,3]);
  t2Sp_a[1]=sum(M[,15]);
  t2Sp_c[1]=sum(1*M[,1]);
  t2Sp_b[1]=sum(1*M[,3]);
  
  t1Sp_a[2]=sum(M[,12]);
  t1Sp_c[2]=sum(1*M[,1]);
  t1Sp_b[2]=sum(1*M[,2]);
  t2Sp_a[2]=sum(M[,12]);
  t2Sp_c[2]=sum(1*M[,1]);
  t2Sp_b[2]=sum(1*M[,2]);
  
  t1Sp_a[3]=sum(M[,18]);
  t1Sp_c[3]=sum(1*M[,1]);
  t1Sp_b[3]=sum(1*M[,4]);
  t2Sp_a[3]=sum(M[,18]);
  t2Sp_c[3]=sum(1*M[,1]);
  t2Sp_b[3]=sum(1*M[,4]);
  
  t1Sp_a[4]=sum(M[,24]);
  t1Sp_c[4]=sum(1*M[,1]);
  t1Sp_b[4]=sum(1*M[,6]);
  t2Sp_a[4]=sum(M[,24]);
  t2Sp_c[4]=sum(1*M[,1]);
  t2Sp_b[4]=sum(1*M[,6]);
  
  t1Sp_a[5]=sum(M[,21]);
  t1Sp_c[5]=sum(1*M[,1]);
  t1Sp_b[5]=sum(1*M[,5]);
  t2Sp_a[5]=sum(M[,21]);
  t2Sp_c[5]=sum(1*M[,1]);
  t2Sp_b[5]=sum(1*M[,5]);
  
  t1Sp_a[6]=sum(MU[,8]);
  t1Sp_c[6]=sum(1*MU[,1]);
  t1Sp_b[6]=sum(1*MU[,2]);
  t2Sp_a[6]=sum(MU[,8]);
  t2Sp_c[6]=sum(1*MU[,1]);
  t2Sp_b[6]=sum(1*MU[,2]);
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(L[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(L[[i]][[1]][L[[i]][[1]][, 2] > 4, 8])
    
   
    for(yy in 1:2)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        flag=1
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
        
        if(0%in%L[[i]][[yy+1]][,3])
        {
          L[[i]][[yy+1]] <- merge_zero_rows(L[[i]][[yy+1]])
        }
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(L[[i]][[yy+1]][,5]==1 & L[[i]][[yy+1]][,9]>0)
        t1blobs_c[yy] = t1blobs_c[yy] + 1*truen_num_blobs
        t1blobs_b[yy] = t1blobs_b[yy] + 1*sum(L[[i]][[yy+1]][,9]>0) 
 
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
    
     
    yy = 2
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
  
   idx = c(1,3,4,2,5,6)
  
  t1blobs = (t1blobs_a / t1blobs_b)[idx]
  t2blobs = (t2blobs_a / t2blobs_c)[idx]
  
  t1CO = (t1CO_a / t1CO_b)[idx]
  t2CO = (t2CO_a / t2CO_c)[idx]
  
  t1Hy = (t1Hy_a / t1Hy_b)[idx]
  t2Hy = (t2Hy_a / t2Hy_c)[idx]
  
  sipo[[kk]]=list( t1blobs,t2blobs,t1CO,t2CO,t1Hy,t2Hy,t1Sp,  t2Sp)
  
  t1Sp = (t1Sp_a / t1Sp_b)[idx]
  t2Sp = (t2Sp_a / t2Sp_c)[idx]
  
  # t1Sp = t1Sp[idx]     
  # t2Sp = t2Sp[idx]  
  
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
  t2Mn = (t2blobs + t2CO + t2Hy + t2Sp) / 4
  
  dblobs = (2 * ( t1blobs) *(  t2blobs))/(( t1blobs)+( t2blobs))
  dCO    =  2 * ( t1CO) *( t2CO)/(( t1CO)+(  t2CO))
  dHy    =  2 * ( t1Hy) *( t2Hy)/(( t1Hy)+( t2Hy))
  dSp    =   2 *( t1Sp) *( t2Sp)/(( t1Sp)+( t2Sp))
  dMn    = (dblobs + dCO + dHy + dSp) / 4
  
 
  
  drfblobs =   (t1blobs_b[idx]+t1blobs_c[idx] -2*t1blobs_a[idx])/(2*t2blobs_c[idx])
  drfCO    =   (t1CO_b[idx]+t1CO_c[idx] -2*t1CO_a[idx])/(2*t2CO_c[idx])
  drfHy    =   (t1Hy_b[idx]+t1Hy_c[idx] -2*t1Hy_a[idx])/(2*t2Hy_c[idx])
  drfSp    =   (t1Sp_b[idx]+t1Sp_c[idx] -2*t1Sp_a[idx])/(2*t2Sp_c[idx])
  drfMn    =   (drfblobs + drfCO + drfHy + drfSp) / 4
 
  
  distanceplot[[kk]] = c(drfblobs,drfCO,drfHy,drfSp,drfMn )
  
  df1 <- data.frame(
      value = c( t1blobs,  t1CO,  t1Hy, t1Sp, t1Mn),
    category = rep(c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
      each = 6)
  )
  
  df1$category <- factor(
    df1$category,
    levels = c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    )
  )
  
  df1$measure <- factor(
    df1$measure,
    levels = c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    )
  )
  
  df2 <- data.frame(
    value = c( t2blobs,  t2CO,  t2Hy, t2Sp, t2Mn),
    category = rep(c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 6)
  )
  
  df2$category <- factor(
    df2$category,
    levels = c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    )
  )
  
  df2$measure <- factor(
    df2$measure,
    levels = c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    )
  )
  
  dfd <- data.frame(
    value = c(dblobs, dCO,dHy,dSp,dMn),
    category = rep(c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 6)
  )
  
  dfd$category <- factor(
    dfd$category,
    levels = c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    )
  )
  
  dfd$measure <- factor(
    dfd$measure,
    levels = c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    )
  )
  
  cols <- c(
    "gray70",
    "#A8DADC", 
    "#274C77",     
    "#A8D5A2",     
    "#4CAF50",    
    "#1B5E20"      
  )
  
  p1<- ggplot(df1,
             aes(x = measure,
                 y = value,
                 fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0("Precision - ",c("Full length","Third of each gene")[kk]),
      x = "",
      y = "",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1))+
    guides(fill = "none")
  
  
  p2<- ggplot(df2,
         aes(x = measure,
             y = value,
             fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0("Recall - ",c("Full length","Third of each gene")[kk]),
      x = "",
      y = "",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1))
  
  
 p3<- ggplot(dfd,
         aes(x = measure,
             y = value,
             fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title =  bquote(F[1]~score~"-"~.(c("Full length", "Third of each gene")[kk])),
      x = "",
      y = "",
      fill = NULL
    ) +
   theme_bw() +
   theme(
     legend.position = "bottom",
     legend.direction = "horizontal",
     legend.box = "horizontal",
     plot.title = element_text(hjust = 0.5)
   ) +
   guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
   coord_cartesian(ylim = c(0, 1))+
   guides(fill = "none")
 
 print(p1+p2+p3)
 
 ggsave(
   paste0(c("Full",'Third')[kk],"data",".pdf"),
   plot = p1+p2+p3,
   width = 11.72,
   height = 5,
   units = "in"
 )
 }
  
  

 
  dfd1 <- data.frame(
    value =   distanceplot[[1]],
    category = rep(c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 6)
  )
  
  dfd1$category <- factor(
    df1$category,
    levels = c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    )
  )
  
  dfd1$measure <- factor(
    dfd1$measure,
    levels = c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    )
  )
  
  
  
  dfd2 <- data.frame(
    value =   distanceplot[[2]],
    category = rep(c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 6)
  )
  
  dfd2$category <- factor(
    dfd2$category,
    levels = c(
      
      "Squirrel",
      "HMS",
      "HWC",
      "L1N",
      "L1N-majority",
      "L1N-unrooted"
      
      
    )
  )
  
  dfd2$measure <- factor(
    df1$measure,
    levels = c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    )
  )
  
  pp1<- ggplot(dfd1,
              aes(x = measure,
                  y = value,
                  fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title =   paste0("(a) Dissimilarity - ", c("Full length", "Third of each gene")[1]),
      x = "",
      y = "",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1))+
    guides(fill = "none")
  
  pp2<- ggplot(dfd2,
               aes(x = measure,
                   y = value,
                   fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title =   paste0("(b) Dissimilarity - ", c("Full length", "Third of each gene")[2]),
      x = "",
      y = "",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(ylim = c(0, 1)) 
  
  cp1<-(pp1 + pp2) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
 
  
  ggsave(
    "DisttoTrue.pdf",
    plot = cp1,
    width = 11.72,
    height = 5,
    units = "in"
  )
  