library(openxlsx)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(MSCquartets)

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2"
location <- file.path(base_dir, "sim_out")
setwd(location)
source(paste0(base_dir,"/auxfiles/Finalcode_bootsupport.R"))
plotsroot <- list()
plotsunroot <- list()

NBfull=matrix(rep(0,2*20),nrow = 20)
NBUfull=matrix(rep(0,2*20),nrow = 20)
for(xx in (1:20)) 
{  
  if(xx<10){xp = paste0('',xx)}else{xp=xx}
  sboot_sampled_lvl1 = read.nets(file=paste0(location,'/rooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[2:101]
  sboot_sampled_tob=list()
  nblobs=0
  n4blobs=0
  
 
  sUboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[2:101]
  sUboot_sampled_tob=list()
  nUblobs=0
  nU4blobs=0
  
  
  for(i in 1:100)
  {
    flag=1
    if(is.null(sboot_sampled_lvl1[[i]]$reticulation)){ flag=0 }
    if(flag==1){ sboot_sampled_tob[[i]]= treeOfBlobs(sboot_sampled_lvl1[[i]])}else{sboot_sampled_tob[[i]]=  sboot_sampled_lvl1[[i]]  }
    deg=ape::degree(sboot_sampled_tob[[i]])
    nblobs = sum( deg[ deg[,1]>3,2]  ) + nblobs
    n4blobs = sum(deg[ deg[,1]>4,2] ) + n4blobs
    
    
    flagU=1
    if(is.null(sUboot_sampled_lvl1[[i]]$reticulation)){ flagU=0 }
    if(flagU==1){ sUboot_sampled_tob[[i]]= treeOfBlobs(sUboot_sampled_lvl1[[i]])}else{sUboot_sampled_tob[[i]]=  sUboot_sampled_lvl1[[i]]  }
    degU=ape::degree(sUboot_sampled_tob[[i]])
    nUblobs = sum( degU[ degU[,1]>3,2] ) + nUblobs
    nU4blobs = sum( degU[ degU[,1]>4,2] ) + nU4blobs
    
  }
  NBfull[xx,1]=nblobs
  NBfull[xx,2]=n4blobs
  
  NBUfull[xx,1]=nUblobs
  NBUfull[xx,2]=nU4blobs
}




NBthird=matrix(rep(0,2*20),nrow = 20)
NBUthird=matrix(rep(0,2*20),nrow = 20)
for(xx in (1:20)) 
{  
  if(xx<10){xp = paste0('',xx)}else{xp=xx}
  sboot_sampled_lvl1 = read.nets(file=paste0(location,'/rooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[2:101]
  sboot_sampled_tob=list()
  nblobs=0
  n4blobs=0
  
  
  sUboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[2:101]
  sUboot_sampled_tob=list()
  nUblobs=0
  nU4blobs=0
  
  
  for(i in 1:100)
  {
    flag=1
    if(is.null(sboot_sampled_lvl1[[i]]$reticulation)){ flag=0 }
    if(flag==1){ sboot_sampled_tob[[i]]= treeOfBlobs(sboot_sampled_lvl1[[i]])}else{sboot_sampled_tob[[i]]=  sboot_sampled_lvl1[[i]]  }
    deg=ape::degree(sboot_sampled_tob[[i]])
    nblobs = sum( deg[ deg[,1]>3,2] ) + nblobs
    n4blobs = sum( deg[ deg[,1]>4,2] ) + n4blobs
    
    flagU=1
    if(is.null(sUboot_sampled_lvl1[[i]]$reticulation)){ flagU=0 }
    if(flagU==1){ sUboot_sampled_tob[[i]]= treeOfBlobs(sUboot_sampled_lvl1[[i]])}else{sUboot_sampled_tob[[i]]=  sUboot_sampled_lvl1[[i]]  }
    degU=ape::degree(sUboot_sampled_tob[[i]])
    nUblobs = sum( degU[ degU[,1]>3,2] ) + nUblobs
    nU4blobs = sum( degU[ degU[,1]>4,2] ) + nU4blobs
    
  }
  NBthird[xx,1]=nblobs
  NBthird[xx,2]=n4blobs
  
  NBUthird[xx,1]=nUblobs
  NBUthird[xx,2]=nU4blobs
}







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

file1=paste0(location,'/EdgeAnalysis_full_nooutgroup.Rdata')
file2=paste0(location,'/EdgeAnalysis_third_nooutgroup.Rdata')
file3=paste0(location,'/EdgeAnalysis_full_outgroup.Rdata')
file4=paste0(location,'/EdgeAnalysis_third_outgroup.Rdata')

load(file1)#M3
load(file2)#M4
load(file3)#M1
load(file4)#M2
####Start
sipo=c()
kk=1
for(kk in 1:2)
{
  NN=list(NBfull,NBthird)[[kk]]
  L =list(result3,result4)[[kk]]
  LU=list(result1,result2)[[kk]]
  
  M=list(MM1,MM2)[[kk]]
  MU=list(MM3,MM4)[[kk]]
  
  
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
 
  t1Sp_a[1]=sum(M[,29]);
  t1Sp_b[1]=sum(1*M[,3]);
 
  
  t1Sp_a[2]=sum(M[,28]);
  t1Sp_b[2]=sum(1*M[,2]);
   
  
  t1Sp_a[3]=sum(M[,30]);
  t1Sp_b[3]=sum(1*M[,4]);
 
  
  t1Sp_a[4]=sum(M[,32]);
  t1Sp_b[4]=sum(1*M[,6]);
  
  
  t1Sp_a[5]=sum(M[,27]);
  t1Sp_b[5]=sum(1*M[,1]);
  
   
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(L[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(L[[i]][[1]][L[[i]][[1]][, 2] > 4, 8])
    
    
    yy=0
    if(L[[i]][[yy+1]][1,1]=='NA'){
      flag=1
    } else {
      L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
      
      
      t1blobs_a[yy+5] = t1blobs_a[yy+5] +   sum(as.numeric(L[[i]][[yy+1]][,4]))*100 #BP vs bootsrap
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
        # bg4 = bg4+1
      }
    }
    
    for(yy in 1:2)
    {
      if(L[[i]][[yy+1]][1,1]=='NA'){
        flag=1
      } else {
        L[[i]][[yy+1]][is.na(L[[i]][[yy+1]])] <- 0
        
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(as.numeric(L[[i]][[yy+1]][,4]))*100 #BP vs bootsrap
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
          # bg4 = bg4+1
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
          t2Hy_a[yy] = t2Hy_a[yy] + truen_num_blobs_ge4*100
          t2Hy_c[yy] = t2Hy_c[yy] + 1*truen_num_blobs_ge4*100
          t1Hy_b[yy] = t1Hy_b[yy] + 1*0
          # bg4 = bg4+1
          
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
          # bg4 = bg4+1
        }
      }
    }
 
  }
  
  idx = c(5,1,3,4,2)
  
  t1blobs = ((t1blobs_b+ rep(sum(NN[,1]),5) - 2*t1blobs_a)/ (2*t1blobs_b))[idx]
  t1CO =((t1CO_b+ rep(sum(NN[,1]),5) - 2*t1CO_a)/ (2*t1CO_b))[idx]
  t1Hy =((t1Hy_b+ rep(sum(NN[,2]),5) - 2*t1Hy_a)/ (2*t1Hy_b))[idx]
  
  t1Sp = (t1Sp_a / (2*t1Sp_b) )[idx]
   
  
  # t1Sp = t1Sp[idx]     
  # t2Sp = t2Sp[idx]  
  
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
 
 
  
   
  
 
  
  df1 <- data.frame(
    value = c( t1blobs,  t1CO,  t1Hy, t1Sp, t1Mn),
    category = rep(c(
      
      "TRUE","Squirrel",
      "HMS",
      "HWC",
      "L1N" 
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 5)
  )
  
  df1$category <- factor(
    df1$category,
    levels = c(
      
      "TRUE","Squirrel",
      "HMS",
      "HWC",
      "L1N" 
      
      
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
   
  cols <- c(
    "gray35",
    "gray70",
    "#A8DADC", 
    "#274C77",     
    "#A8D5A2"     
  )
  
  p1<- ggplot(df1,
              aes(x = measure,
                  y = value,
                  fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0(c("(a) ","(c) ")[kk], "Dissimilarity to Input networks - Rooted - ",c("Full length","Third of each gene")[kk]),
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
    coord_cartesian(ylim = c(0,1.2)) 
  plotsroot[[kk]] <- p1
}











########### ########### ########### 
###########   UNROOTED  ########### 
########### ########### ########### 


sipo=c()
kk=1
for(kk in 1:2)
{
  NN=list(NBUfull,NBUthird)[[kk]]
  
  LU=list(result1,result2)[[kk]]
  
  
  MU=list(MM3,MM4)[[kk]]
  
  
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
  
  t1Sp_a[1]=sum(MU[,19]);
  t1Sp_b[1]=sum(1*MU[,1]);
  
  
  t1Sp_a[2]=sum(MU[,18]);
  t1Sp_b[2]=sum(1*MU[,2]);
  
  
  t1Sp_a[3]=sum(MU[,17]);
  t1Sp_b[3]=sum(1*MU[,4]);
  
   
  
  
  for(i in 1:20)
  {
    truen_num_blobs      = length(LU[[i]][[1]][,1])
    truen_num_blobs_ge4  = length(LU[[i]][[1]][LU[[i]][[1]][, 2] > 4, 8])
    
    
    yy=0
    if(LU[[i]][[yy+1]][1,1]=='NA'){
      flag=1
    } else {
      LU[[i]][[yy+1]][is.na(LU[[i]][[yy+1]])] <- 0
      
      
      t1blobs_a[yy+3] = t1blobs_a[yy+3] +   sum(as.numeric(LU[[i]][[yy+1]][,4]))*100 #BP vs bootsrap
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
        # bg4 = bg4+1
      }
    }
    
    for(yy in 1:2)
    {
      if(LU[[i]][[yy+1]][1,1]=='NA'){
        flag=1
      } else {
        LU[[i]][[yy+1]][is.na(LU[[i]][[yy+1]])] <- 0
        
        
        t1blobs_a[yy] = t1blobs_a[yy] +   sum(as.numeric(LU[[i]][[yy+1]][,4]))*100 #BP vs bootsrap
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
          # bg4 = bg4+1
        }
      }
    }
 
  }
  
  idx = c(3,1,2)
  
  t1blobs = ((t1blobs_b+ rep(sum(NN[,1]),3) - 2*t1blobs_a)/ (2*t1blobs_b))[idx]
  
  
  t1CO =((t1CO_b+ rep(sum(NN[,1]),3) - 2*t1CO_a)/ (2*t1CO_b))[idx]
  
  
  t1Hy =((t1Hy_b+ rep(sum(NN[,2]),3) - 2*t1Hy_a)/ (2*t1Hy_b))[idx]
  
  
  t1Sp = (t1Sp_a / (2*t1Sp_b) )[idx]
  
  
  # t1Sp = t1Sp[idx]     
  # t2Sp = t2Sp[idx]  
  
  t1Mn = (t1blobs + t1CO + t1Hy + t1Sp) / 4
  
  
  
  
  
  
  
  df1 <- data.frame(
    value = c( t1blobs,  t1CO,  t1Hy, t1Sp, t1Mn),
    category = rep(c(
      
      "TRUE","Squirrel", "L1N" 
      
      
    ), 5),
    measure = rep(c(
      "Blobs",
      "Cycles",
      "Hybrids B>4",
      "Splits",
      "Mean"
    ),
    each = 3)
  )
  
  df1$category <- factor(
    df1$category,
    levels = c(
      
      "TRUE","Squirrel",
       
      "L1N" 
      
      
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
  
  cols <- c(
    "gray35",
    "gray70",
     
    "#A8D5A2"     
  )
  
  p1<- ggplot(df1,
              aes(x = measure,
                  y = value,
                  fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste0(c("(b) ","(d) ")[kk], "Dissimilarity to Input networks - Unrooted - ",c("Full length","Third of each gene")[kk]),
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
    coord_cartesian(ylim = c(0,1.2)) 
  plotsunroot[[kk]] <- p1
}

# grid.arrange(
#   grobs = c(rbind(plotsroot, plotsunroot)),
#   ncol = 2
# )
 
grid.arrange(
  grobs = c(cbind(plotsroot[[1]], plotsunroot[[1]])),
  ncol = 2
)

combined_plot <- grid.arrange(
  grobs = c(plotsroot[[1]], plotsunroot[[1]]),
  ncol = 2
)

ggsave(
  "DistFull.pdf",
  plot = combined_plot,
  width = 11.72,
  height = 5,
  units = "in"
)

combined_plot <- grid.arrange(
  grobs = c(cbind(plotsroot[[2]], plotsunroot[[2]])),
  ncol = 2
)

ggsave(
  "DistHalf.pdf",
  plot = combined_plot,
  width = 11.72,
  height = 5,
  units = "in"
)
