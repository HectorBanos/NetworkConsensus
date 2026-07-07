library(stringr) 
library(ape)
library(SiPhyNetwork)
library(MSCquartets)
base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2"
location <- file.path(base_dir, "sim_out")
source(paste0(base_dir,"/auxfiles/Finalcode_bootsupport.R"))

TP_bp =c()
TP_ts =c()
FP_bp =c()
FP_ts =c()
dist=c()
dist2=c()
setwd(paste0(base_dir,"/auxfiles"))
networks = readLines("networksjl_outgroup_adjusted_labelledv2.txt")
 
tobs=  readLines('lvl1networks_tob_letters.txt')
things=list()
for(xx in 1:20)
{

         things[[xx]]=list()
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

         at = blob_dist_boot(sToB,sboot_sampled_tob)
         st = blob_dist_boot(sqToB,sboot_sampled_tob)
         ct = blob_dist_boot(TrueToB,sboot_sampled_tob)

         ff=blob_dist_boot(sToB,list(TrueToB,TrueToB))
         sff=blob_dist_boot(sqToB,list(TrueToB,TrueToB))

         ac = cycle_dist_boot(sLvl1,sboot_sampled_lvl1)
         sc = cycle_dist_boot(sqLvl1,sboot_sampled_lvl1)
         cc = cycle_dist_boot(TrueNet,sboot_sampled_lvl1)

         dd=cycle_dist_boot(sLvl1,list(TrueNet,TrueNet))
         sdd=cycle_dist_boot(sqLvl1,list(TrueNet,TrueNet))


         ah = hybrid_dist_boot(sLvl1,sboot_sampled_lvl1)
         sh = hybrid_dist_boot(sqLvl1,sboot_sampled_lvl1)
         ch = hybrid_dist_boot(TrueNet,sboot_sampled_lvl1)


         hh = hybrid_dist_boot(sLvl1,list(TrueNet,TrueNet))
         shh = hybrid_dist_boot(sqLvl1,list(TrueNet,TrueNet))



         things[[xx]][[1]]=at
         things[[xx]][[2]]=st
         things[[xx]][[3]]=ct
         things[[xx]][[4]]=ff
         things[[xx]][[5]]=sff
         things[[xx]][[6]]=ac
         things[[xx]][[7]]=sc
         things[[xx]][[8]]=cc
         things[[xx]][[9]]=dd
         things[[xx]][[10]]=sdd
         things[[xx]][[11]]=ah
         things[[xx]][[12]]=sh
         things[[xx]][[13]]=ch
         things[[xx]][[14]]=hh
         things[[xx]][[15]]=shh


         print(xx)

}

  result=list()
  for(xx in 1:20)
  {
    result[[xx]]=list()
    result[[xx]][[1]] <- merge_three_trees(things[[xx]][[3]], things[[xx]][[8]],things[[xx]][[13]])
    # View(result[,-1])

    result[[xx]][[2]] <- merge_six_trees(things[[xx]][[2]], things[[xx]][[5]],things[[xx]][[7]], things[[xx]][[10]],things[[xx]][[12]],things[[xx]][[15]])
    # View(result[,-1])

    result[[xx]][[3]] <- merge_six_trees(things[[xx]][[1]], things[[xx]][[4]],things[[xx]][[6]], things[[xx]][[9]],things[[xx]][[11]],things[[xx]][[14]])
    # View(result[,-1])


  }
  save(result,file=paste0(location,'/BlobAnalysis_full_nooutgroup.Rdata'))



  
  things=list()
  for(xx in 1:20)
  {      
    
    things[[xx]]=list()
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
    
    at = blob_dist_boot(sToB,sboot_sampled_tob)
    st = blob_dist_boot(sqToB,sboot_sampled_tob)
    ct = blob_dist_boot(TrueToB,sboot_sampled_tob)
    
    ff=blob_dist_boot(sToB,list(TrueToB,TrueToB))
    sff=blob_dist_boot(sqToB,list(TrueToB,TrueToB))
    
    ac = cycle_dist_boot(sLvl1,sboot_sampled_lvl1)
    sc = cycle_dist_boot(sqLvl1,sboot_sampled_lvl1)
    cc = cycle_dist_boot(TrueNet,sboot_sampled_lvl1)
    
    dd=cycle_dist_boot(sLvl1,list(TrueNet,TrueNet))
    sdd=cycle_dist_boot(sqLvl1,list(TrueNet,TrueNet))
    
    
    ah = hybrid_dist_boot(sLvl1,sboot_sampled_lvl1)
    sh = hybrid_dist_boot(sqLvl1,sboot_sampled_lvl1)
    ch = hybrid_dist_boot(TrueNet,sboot_sampled_lvl1)
    
    
    hh = hybrid_dist_boot(sLvl1,list(TrueNet,TrueNet))
    shh = hybrid_dist_boot(sqLvl1,list(TrueNet,TrueNet))
    
    
    
    things[[xx]][[1]]=at
    things[[xx]][[2]]=st
    things[[xx]][[3]]=ct
    things[[xx]][[4]]=ff
    things[[xx]][[5]]=sff
    things[[xx]][[6]]=ac
    things[[xx]][[7]]=sc
    things[[xx]][[8]]=cc
    things[[xx]][[9]]=dd
    things[[xx]][[10]]=sdd
    things[[xx]][[11]]=ah
    things[[xx]][[12]]=sh
    things[[xx]][[13]]=ch
    things[[xx]][[14]]=hh
    things[[xx]][[15]]=shh
    
    print(xx)
    
  }      
 
  result=list()
  for(xx in 1:20)
  {
    result[[xx]]=list()
  result[[xx]][[1]] <- merge_three_trees(things[[xx]][[3]], things[[xx]][[8]],things[[xx]][[13]])
    # View(result[,-1])

  result[[xx]][[2]] <- merge_six_trees(things[[xx]][[2]], things[[xx]][[5]],things[[xx]][[7]], things[[xx]][[10]],things[[xx]][[12]],things[[xx]][[15]])
    # View(result[,-1])

  result[[xx]][[3]] <- merge_six_trees(things[[xx]][[1]], things[[xx]][[4]],things[[xx]][[6]], things[[xx]][[9]],things[[xx]][[11]],things[[xx]][[14]])
   # View(result[,-1])
 
  }
  save(result,file=paste0(location,'/BlobAnalysis_third_nooutgroup.Rdata'))
  
