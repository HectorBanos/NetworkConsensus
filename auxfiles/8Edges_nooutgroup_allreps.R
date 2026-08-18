#########################################################################
#################### Process edge info unrooted case ####################
#########################################################################

### This script was written by a human and optimized using Claude LLM

library(ape)
library(SiPhyNetwork)
library(grid)
library(gridExtra)
library(TreeDist)
library(MSCquartets)
library('Quartet')

 
base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"

source(paste0(base_dir,"/auxfiles/Finalcode_bootsupport.R"))

edge_dist <- function(tree, boot_tree) {
  M1 <- split_matrix(tree)[-1,]
  names1  <- sapply(1:nrow(M1), function(i) {
    row <- M1[i, ]
    target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })
  M2 <- split_matrix(boot_tree)[-1,]
  names2  <- sapply(1:nrow(M2), function(i) {
    row <- M2[i, ]
    target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })

  matches <- intersect(names1, names2)
  a=length(matches)
  only_a <- setdiff(names1, names2)
  b=length(only_a)
  only_b <- setdiff(names2, names1)
  c= length(only_b)

  edge_dist <-c(a,b,c)
}

edge_supp <- function(tree, boot_trees) {

  M1 <- split_matrix(tree)[-1,]
  names1  <- sapply(1:nrow(M1), function(i) {
    row <- M1[i, ]
    target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })
  edges=rep(0,length(names1))
  for(j in 1:length(names1))
  {
    nam=names1[j]
    for(i in 1:length(boot_trees))
    {
      M2 <- split_matrix(boot_trees[[i]])[-1,]
      names2  <- sapply(1:nrow(M2), function(i) {
        row <- M2[i, ]
        target <- row[1]
        paste(names(row[row == target]), collapse = ",")
      })
      edges[j]= nam%in% names2+edges[j]
    }

  }

  edge_supp <-   names1[which(edges>50)]

}

RF.DIFM<- function(boot_tree ,tree, boot_trees)
{
  names2=edge_supp(tree, boot_trees)
  M1 <- split_matrix(boot_tree)[-1,]
  names1  <- sapply(1:nrow(M1), function(i) {
    row <- M1[i, ]
    target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })

  matches <- intersect(names1, names2)
  a=length(matches)
  only_a <- setdiff(names1, names2)
  b=length(only_a)
  only_b <- setdiff(names2, names1)
  c= length(only_b)

  RF.DIFM<-c( length(names1)+length(names2) -2*a , a,b,c)

}

setwd(paste0(base_dir,"/auxfiles"))

networks = readLines("networksjl_outgroup_adjusted_labelled.txt")
loopsize = 20
tobs=  readLines('lvl1networks_tob_letters.txt')

n_reps <- 10

for (rep_i in 1:n_reps) {

  message("=== Edges nooutgroup: replicate ", rep_i, " ===")

  ## ============= FULL block (unrootedfull<rep_i>) =============
  location <- file.path(base_dir, "sim_out", paste0("unrootedfull", rep_i))
  setwd(location)

  MM=matrix(rep(0,20*20),nrow = 20)
  for(xx in (1:loopsize) )
  {

    TrueNet= read.evonet(text=networks[xx])
    newick <- networks[xx]
    TrueToB=treeOfBlobs( TrueNet)
    TrueToB$edge.length = rep(1,length(TrueToB$edge.length))
    if(xx<10){xp = paste0('',xx)}else{xp=xx}

    sLvl1=read.evonets(file=paste0(location,'/Consensus_nooutgroup_full.txt'))[[xx]]
    sqLvl1=read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[1][[1]]
    sToB=treeOfBlobs(sLvl1)
    sqToB=treeOfBlobs(sqLvl1)
    sToB$edge.length = rep(1,length(sToB$edge.length))

    sboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[2:101]
    sboot_sampled_tob=list()
    for(i in 1:100)
    {
      flag=1
      if(is.null(sboot_sampled_lvl1[[i]]$reticulation)){ flag=0 }
      if(flag==1){ sboot_sampled_tob[[i]]= treeOfBlobs(sboot_sampled_lvl1[[i]])}else{sboot_sampled_tob[[i]]=  sboot_sampled_lvl1[[i]]  }
    }

    srf=RF.dist(TrueToB,sToB)
    sqrf=RF.dist(TrueToB,sqToB)
    mrf= RF.DIFM(TrueToB, sToB,sboot_sampled_tob)[1]

    numbran=dim(TrueToB$edge)[1] -21
    numbrans=dim(sToB$edge)[1] -21
    numbransq=dim(sqToB$edge)[1] -21
    numbranm=length(edge_supp(sToB,sboot_sampled_tob))

    sed=edge_dist(TrueToB,sToB)
    sqed=edge_dist(TrueToB,sqToB)
    med=RF.DIFM(TrueToB, sToB,sboot_sampled_tob)[2:4]

    srfbs=0
    sqrfbs=0
    hrfbs=0
    mrfbs=0
    Trfbs= 0

    for(i in 1:100)
    {
      srfbs=RF.dist(sboot_sampled_tob[[i]],sToB)/100+srfbs
      sqrfbs=RF.dist(sboot_sampled_tob[[i]],sqToB)/100+sqrfbs
      mrfbs= RF.DIFM(sboot_sampled_tob[[i]], sToB,sboot_sampled_tob)[1]/100+mrfbs
      Trfbs=RF.dist(sboot_sampled_tob[[i]],TrueToB)/100+Trfbs
    }

    MM[xx,]=c(numbran,numbrans,numbransq,numbranm,srf,sqrf,mrf,sed,sqed,med,Trfbs,srfbs,sqrfbs,mrfbs)
  }

  MM3=MM
  save(MM3,file=paste0(location,'/EdgeAnalysis_full_nooutgroup.Rdata'))

  ## ============= THIRD block (unrootedthird<rep_i>) =============
  location <- file.path(base_dir, "sim_out", paste0("unrootedthird", rep_i))
  setwd(location)

  MM=matrix(rep(0,20*20),nrow = 20)
  for(xx in (1:loopsize) )
  {

    TrueNet= read.evonet(text=networks[xx])
    newick <- networks[xx]
    TrueToB=treeOfBlobs( TrueNet)
    TrueToB$edge.length = rep(1,length(TrueToB$edge.length))
    if(xx<10){xp = paste0('',xx)}else{xp=xx}

    flags=1
    sLvl1=readLines( paste0(location,'/Consensus_nooutgroup_third.txt'))[xx]
    if(str_detect(sLvl1, "#"))
    {  sLvl1=read.evonet(text=sLvl1)
    }else{ sLvl1=read.tree(text=sLvl1)
    flags=0;
    }

    sqLvl1=read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[1][[1]]
    if(flags==1){ sToB=treeOfBlobs(sLvl1)}else{sToB= unroot(sLvl1) }
    sqToB=treeOfBlobs(sqLvl1)
    sToB$edge.length = rep(1,length(sToB$edge.length))

    sboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[2:101]
    sboot_sampled_tob=list()
    for(i in 1:100)
    {
      flag=1
      if(is.null(sboot_sampled_lvl1[[i]]$reticulation)){ flag=0 }
      if(flag==1){ sboot_sampled_tob[[i]]= treeOfBlobs(sboot_sampled_lvl1[[i]])}else{sboot_sampled_tob[[i]]=  sboot_sampled_lvl1[[i]]  }
    }

    srf=RF.dist(TrueToB,sToB)
    sqrf=RF.dist(TrueToB,sqToB)
    mrf= RF.DIFM(TrueToB, sToB,sboot_sampled_tob)[1]

    numbran=dim(TrueToB$edge)[1] -21
    numbrans=dim(sToB$edge)[1] -21
    numbransq=dim(sqToB$edge)[1] -21
    numbranm=length(edge_supp(sToB,sboot_sampled_tob))

    sed=edge_dist(TrueToB,sToB)
    sqed=edge_dist(TrueToB,sqToB)
    med=RF.DIFM(TrueToB, sToB,sboot_sampled_tob)[2:4]

    srfbs=0
    sqrfbs=0
    hrfbs=0
    mrfbs=0
    Trfbs= 0

    for(i in 1:100)
    {
      srfbs=RF.dist(sboot_sampled_tob[[i]],sToB)/100+srfbs
      sqrfbs=RF.dist(sboot_sampled_tob[[i]],sqToB)/100+sqrfbs
      mrfbs= RF.DIFM(sboot_sampled_tob[[i]], sToB,sboot_sampled_tob)[1]/100+mrfbs
      Trfbs=RF.dist(sboot_sampled_tob[[i]],TrueToB)/100+Trfbs
    }

    MM[xx,]=c(numbran,numbrans,numbransq,numbranm,srf,sqrf,mrf,sed,sqed,med,Trfbs,srfbs,sqrfbs,mrfbs)
  }

  MM4=MM
  save(MM4,file=paste0(location,'/EdgeAnalysis_third_nooutgroup.Rdata'))

  message("Completed replicate ", rep_i, " of ", n_reps)
}

message("All ", n_reps, " replicates complete (edges, nooutgroup).")
