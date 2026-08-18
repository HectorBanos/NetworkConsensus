xnblobs=0
xn4blobs =0

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2"
location <- file.path(base_dir, "sim_out")
setwd(location)
for(xx in (1:20)) 
{  
  if(xx<10){xp = paste0('',xx)}else{xp=xx}
  sboot_sampled_lvl1 = read.nets(file=paste0(location,'/rooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[1]
  suboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'.txt'))[1]
  stboot_sampled_lvl1 = read.nets(file=paste0(location,'/rooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[1]
  stuboot_sampled_lvl1 = read.nets(file=paste0(location,'/unrooted_network_newick_concatalignRealisticSim_400gts_Net',xp,'_third.txt'))[1]
  # 
  i=1
    flag=1
    if(is.null(sboot_sampled_lvl1[[i]]$reticulation)){ flag=0 }
    if(flag==1){ sboot_sampled_tob[[i]]= treeOfBlobs(sboot_sampled_lvl1[[i]])}else{sboot_sampled_tob[[i]]=  sboot_sampled_lvl1[[i]]  }
    deg=ape::degree(sboot_sampled_tob[[i]])
    xnblobs = sum( deg[ deg[,1]>3,2]  ) + xnblobs
    xn4blobs = sum( deg[ deg[,1]>4,2]  ) + xn4blobs
   
    
}
xnblobs/20
xn4blobs/20


# Rooted full
# > t1Mn
# [1] 0.2692223 0.2109775 0.3082724 0.3753080 0.1764533
# > 0.2692223-0.2109775
# [1] 0.0582448

# #Rooted third
# > t1Mn
# [1] 0.3214072 0.2745655 0.5511568 0.4403982 0.2316099
# > 0.3214072-0.2745655
# [1] 0.0468417

#Unrooted full
# >  t1Mn
# [1] 0.2843682 0.2309212 0.1869746
# > 0.2843682-0.2309212
# [1] 0.053447

# Unrooted thrd
# >  t1Mn
# [1] 0.3302704 0.2888582 0.2421767
# > 0.3302704-0.2888582
# [1] 0.0414122
