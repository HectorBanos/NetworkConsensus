##### Code to compute the BP and TS distance for Blob, Cricular, and Hybrid support
library(RcppHungarian)
library(igraph)
library(ape)
library(SiPhyNetwork)
library(phangorn)
library(tibble)
library(stringr)
 


### Read a file with a collectiopn of networks *like read.tree* but for networks
#Input: file = the name of a file with networks in Newick format, where each line is a network


#############################
#############################
#####  Hybrid Distance  #####
#############################
#############################
#############################
#############################
### Compute the Hybrid TS and BP distance for a pair of networks (with the first one being the reference and teh second the comparative network). For all the blobs in the reference network, ut computes the hybrid distance to all the blobs in the comparative ntwk and outputs the minimum
#Input: net1=a network in evonet format,  net2=a network in evonet format
 hybrid_dist <- function(net1, net2,min.length=4) {
  circord <- cycle_circord(net1) 
  circordst <-circord 

  #### Need to fix the 4 cycle
  h1 <- net1$reticulation
  h2 <- net2$reticulation
  ch1 <- lapply(h1[, 2], function(node) {
    net1$tip.label[unlist(Descendants(net1, node, type = "tips"))]
  })
  
  
  problematic <- get_taxa_below_2cycles(net1)$taxa
  if(length(problematic)>0)
  {
  ch1<- ch1[ ch1 !=  problematic]
  }
  circord <- cycle_circordv2(net1)
  circord  <- circord [order(sapply(circord , length))]
  circordsmall  <-  circord[unlist(lapply(circord,length))<4]
  tinymisfits<- matching_hybrid_and_cycles(  circordsmall ,ch1)
  ch1<-ch1[!tinymisfits>0]
  circord <- cycle_circord(net1)
  hybrid_order <- matching_hybrid_and_cycles(circord,ch1)
  if(sum(hybrid_order==0)>0)
  {
    ch1 <- ch1[- which(hybrid_order==0)]
  }
  circord<-circord[hybrid_order]
  labelorder <- cycle_circord_label(circordst)[hybrid_order]
  # ch1 <- lapply(h1[, 2], function(node) {
  #   net1$tip.label[unlist(Descendants(net1, node, type = "tips"))]
  # })
  
  # hybrid_order <- matching_hybrid_and_cycles(circord,ch1)
  # 
  # circord <- circord[hybrid_order]
  if(!is.null(h2))
  {
  # Precompute descendant tip sets for each reticulation

  
  ch2 <- lapply(h2[, 2], function(node) {
    net2$tip.label[unlist(Descendants(net2, node, type = "tips"))]
  })
  
  problematic <- get_taxa_below_2cycles(net2)$taxa
  if(length(problematic)>0)
  {
    ch2<- ch2[ ch2 !=  problematic]
  }
  
  # Vectorized score calculation
  score <- vapply(seq_along(ch1), function(i) {
    len1 <- length(ch1[[i]])
    
    dists <- vapply(ch2, function(c2) {
      len2 <- length(c2)
      inter <- length(intersect(ch1[[i]], c2))
      (len1 + len2 - 2 * inter) / (len1 + len2)
    }, numeric(1))
    
    min(dists)
  }, numeric(1)) }else{score = 1}
  return(cbind(labelorder,score))
 }
 
 
 hybrid_distv2 <- function(net1, net2,min.length=4) {
   circord <- cycle_circordv2(net1) 
   circordst <-circord 
   
   #### Need to fix the 4 cycle
   h1 <- net1$reticulation
   h2 <- net2$reticulation
   ch1 <- lapply(h1[, 2], function(node) {
     net1$tip.label[unlist(Descendants(net1, node, type = "tips"))]
   })
   
   
   problematic <- get_taxa_below_2cycles(net1)$taxa
   if(length(problematic)>0)
   {
     ch1<- ch1[ ch1 !=  problematic]
   }
   circord <- cycle_circordv2(net1)
   circord  <- circord [order(sapply(circord , length))]
   circord <- cycle_circordv2(net1)
   hybrid_order <- matching_hybrid_and_cycles(circord,ch1)
    
   circord<-circord[hybrid_order]
   labelorder <- cycle_circord_label(circordst)[hybrid_order]
    
   if(!is.null(h2))
   {
     # Precompute descendant tip sets for each reticulation
     
     
     ch2 <- lapply(h2[, 2], function(node) {
       net2$tip.label[unlist(Descendants(net2, node, type = "tips"))]
     })
     
     problematic <- get_taxa_below_2cycles(net2)$taxa
     if(length(problematic)>0)
     {
       ch2<- ch2[ ch2 !=  problematic]
     }
     
     # Vectorized score calculation
     score <- vapply(seq_along(ch1), function(i) {
       len1 <- length(ch1[[i]])
       
       dists <- vapply(ch2, function(c2) {
         len2 <- length(c2)
         inter <- length(intersect(ch1[[i]], c2))
         (len1 + len2 - 2 * inter) / (len1 + len2)
       }, numeric(1))
       
       min(dists)
     }, numeric(1)) }else{score = 1}
   return(cbind(labelorder,score))
 }
 
 
 #Computes the bootstrap support 
 hybrid_dist_boot <- function(lvl1, bootlvl1) {
   
    problematic <- get_taxa_below_2cycles(lvl1)$taxa
 
   
   h <- lvl1$reticulation 
   
  if(!is.null(h))
  {
   ch1 <- lapply(h[, 2], function(node) {
     lvl1$tip.label[unlist(Descendants(lvl1, node, type = "tips"))]
   })
   if(length(problematic)>0){
   ch1<- ch1[ ch1 !=  problematic]
   }
 
   
   # Compute all scores with lapply (vectorized + cleaner)
   scores <- do.call(rbind, lapply(bootlvl1, function(net2) hybrid_dist(lvl1, net2)))
   
   circord <- cycle_circordv2(lvl1)
   circord  <- circord [order(sapply(circord , length))]
   
   circordsmall  <-  circord[unlist(lapply(circord,length))<4]
   tinymisfits<- matching_hybrid_and_cycles(  circordsmall ,ch1)
   ch1<-ch1[!tinymisfits>0]
   
   circord  <-  circord[unlist(lapply(circord,length))>3]
   hybrid_order <- matching_hybrid_and_cycles(circord,ch1)
   
   
   
   if(sum(hybrid_order==0)>0)
   {
     ch1 <- ch1[- which(hybrid_order==0)]
   }
   circord<-circord[hybrid_order]
   inull <-  which(sapply(circord, is.null)==T)
   circord <- circord[!sapply(circord, is.null)]
   i3cyc <-  which(sapply(circord, length)==3)
   circord <- circord[!(sapply(circord, length)==3)]
   labelorder <- cycle_circord_label(circord)
   if(length(inull)>0)
   {ch1 <- ch1[-inull];}
   if(length(i3cyc)>0)
   {  ch1 <- ch1[- i3cyc]}
    
   # honrd <- vapply(ch1, function(x) {
   #   match(TRUE, vapply(circord, function(y) all( x == y ), logical(1)))
   # }, integer(1))
   
   honrd=c()
   for(i in 1:length(ch1))
   {
     for(j in 1:length(circord))
     {
       for(k in 1:length(circord[[j]]))
        if( setequal(ch1[[i]] ,circord[[j]][[k]])  )
        {
          honrd[i]=j
        }
     }
   }
 
   
   
   hybridscores <- c("Circular order","Cycle size","Hybrids","TS","BP")
   for(i in 1:length(labelorder))
       {
         score = as.numeric(scores[scores[,1]==labelorder[i],][,2])
         
         hybridscores=rbind(hybridscores, c(labelorder[i],length(circord[[i]]),paste0(ch1[[honrd[i]]],collapse = '*') ,round(1 - mean(score),3),round((sum(score==0))/length(score),3 )) )
   }
   
   hybridscores <- hybridscores[-1,]
   if(is.vector(hybridscores)){ hybridscores<- matrix(hybridscores, ncol = 5)}
   if(is.matrix(hybridscores)){
   hybridscores <- `colnames<-`(hybridscores,c("Circular order","Cycle size","Hybrids","TS","BP")) 
   }else{
     names(hybridscores) <- c("Circular order","Cycle size","Hybrids","TS","BP")
   }
   
     # hybridscores <- hybridscores[  (hybridscores[,2]=="3") + (hybridscores[,2]=="4")<1 ,  ]
     hybridscores <- hybridscores[  (hybridscores[,2]=="3")  <1 ,  ]
     
     
   if( is.null(dim(hybridscores)) && length(hybridscores)>0  ){hybridscores <- t(hybridscores)
   
   }else if(!is.null(dim(hybridscores))){hybridscores<-hybridscores}else{hybridscores=('Input network has no hybrids of length >4')}
     
   return(hybridscores)
  }else{return('Input network has no hybrids')}
   # Bootstrap consensus scores
 }
 
 
 
 hybrid_dist_bootv2 <- function(lvl1, bootlvl1) {
   
   problematic <- get_taxa_below_2cycles(lvl1)$taxa
   
   
   h <- lvl1$reticulation 
   
   if(!is.null(h))
   {
     ch1 <- lapply(h[, 2], function(node) {
       lvl1$tip.label[unlist(Descendants(lvl1, node, type = "tips"))]
     })
     if(length(problematic)>0){
       ch1<- ch1[ ch1 !=  problematic]
     }
     
     
     # Compute all scores with lapply (vectorized + cleaner)
     scores <- do.call(rbind, lapply(bootlvl1, function(net2) hybrid_distv2(lvl1, net2,min.length = 3)))
     
     circord <- cycle_circordv2(lvl1)
     circord  <- circord [order(sapply(circord , length))]
     
     # circordsmall  <-  circord[unlist(lapply(circord,length))<4]
     # tinymisfits<- matching_hybrid_and_cycles(  circordsmall ,ch1)
     # ch1<-ch1[!tinymisfits>0]
     
     # circord  <-  circord[unlist(lapply(circord,length)) ]
     hybrid_order <- matching_hybrid_and_cycles(circord,ch1)
     
     
     
     # if(sum(hybrid_order==0)>0)
     # {
     #   ch1 <- ch1[- which(hybrid_order==0)]
     # }
     # circord<-circord[hybrid_order]
     # inull <-  which(sapply(circord, is.null)==T)
     # circord <- circord[!sapply(circord, is.null)]
     # i3cyc <-  which(sapply(circord, length)==3)
     # circord <- circord[!(sapply(circord, length)==3)]
       labelorder <- cycle_circord_label(circord)
     # if(length(inull)>0)
     # {ch1 <- ch1[-inull];}
     # if(length(i3cyc)>0)
     # {  ch1 <- ch1[- i3cyc]}
     
     # honrd <- vapply(ch1, function(x) {
     #   match(TRUE, vapply(circord, function(y) all( x == y ), logical(1)))
     # }, integer(1))
     
     honrd=c()
     for(i in 1:length(ch1))
     {
       for(j in 1:length(circord))
       {
         for(k in 1:length(circord[[j]]))
           if( setequal(ch1[[i]] ,circord[[j]][[k]])  )
           {
             honrd[i]=j
           }
       }
     }
     
     
     
     hybridscores <- c("Circular order","Cycle size","Hybrids","TS","BP")
     for(i in 1:length(labelorder))
     {
       score = as.numeric(scores[scores[,1]==labelorder[i],][,2])
       
       hybridscores=rbind(hybridscores, c(labelorder[i],length(circord[[i]]),paste0(ch1[[honrd[i]]],collapse = '*') ,round(1 - mean(score),3),round((sum(score==0))/length(score),3 )) )
     }
     
     hybridscores <- hybridscores[-1,]
     if(is.vector(hybridscores)){ hybridscores<- matrix(hybridscores, ncol = 5)}
     if(is.matrix(hybridscores)){
       hybridscores <- `colnames<-`(hybridscores,c("Circular order","Cycle size","Hybrids","TS","BP")) 
     }else{
       names(hybridscores) <- c("Circular order","Cycle size","Hybrids","TS","BP")
     }
     
     # hybridscores <- hybridscores[  (hybridscores[,2]=="3") + (hybridscores[,2]=="4")<1 ,  ]
      # hybridscores <- hybridscores[  (hybridscores[,2]=="3")  <1 ,  ]
     
     
     if( is.null(dim(hybridscores)) && length(hybridscores)>0  ){hybridscores <- t(hybridscores)
     
     }else if(!is.null(dim(hybridscores))){hybridscores<-hybridscores}else{hybridscores=('Input network has no hybrids of length >4')}
     
     return(hybridscores)
   }else{return('Input network has no hybrids')}
   # Bootstrap consensus scores
 }
 
 cycle_circord <- function(lvl1) {
   taxa <- lvl1$tip.label
   il1 <- igraph::as.igraph(lvl1)
   
  
   cut_edges <- igraph::bridges(il1)
   non_cut_edges <- base::setdiff(E(il1), cut_edges)
   
    
   indeg <- igraph::degree(il1, mode = "in")
   nodes_indegree_2 <- V(il1)[indeg == 2]
 
   g_sub <- igraph::as_undirected(
     igraph::subgraph_from_edges(il1, eids = non_cut_edges, delete.vertices = TRUE),
     mode = "collapse"
   )
   comp <- igraph::components(g_sub)$membership
   
   
   ord <- lapply(seq_along(nodes_indegree_2), function(i) {
     cycs <- names(comp[comp == i])
     cycle <- igraph::induced_subgraph(g_sub, cycs)
     
     girth_info <- igraph::girth(cycle, circle = TRUE)
     if (girth_info$girth == Inf) return(NULL)
     
     lapply(seq_along(girth_info$circle), function(j) {
       g3 <- igraph::delete_vertices(il1, names(girth_info$circle[j]))
       comp3 <- igraph::components(g3)$membership
       
       
       index <- base::setdiff(c(1, 2), unique(comp3[names(girth_info$circle)[-j]]))
       intersect(names(comp3[comp3 == index]), taxa)
     })
   })
   
   return(ord)
 }
 
 cycle_circord_label <- function(ord) {
   ord_label <- vapply(ord, function(group) {
     paste0(vapply(group, function(x) paste0(x, collapse = "*"), character(1L)), collapse = " / ")
   }, character(1L))
   return(ord_label)
 }
 
 matching_hybrid_and_cycles <- function(ord,hybrid)
{
 
   hyb_ord_map <- vapply(seq_along(hybrid), function(j) {
     match_idx <- which(vapply(ord, function(cycle) {
       any(vapply(cycle, identical, logical(1), y = hybrid[[j]]))
     }, logical(1)))
     if (length(match_idx) > 0) match_idx[1] else 0L
   }, integer(1))
  
   return(hyb_ord_map)
 }
  
 
   hybrid_dist_boot_boot  <- function(bootlvl1,lvl1,lim ) {
   
   pic=list()
   for(i in 1:length(bootlvl1))
   {
     pic[[i]]=hybrid_dist_boot(bootlvl1[[i]],bootlvl1)
   }
   combined <- do.call(rbind,pic)
   
   # Step 2: Remove duplicate rows based on the first column
   unique_rows <- combined[!duplicated(combined[, 1]), ]
   
   summ<- cycle_circord(lvl1)
   summ<-cycle_circord_label(summ)
   unique_rows <-   unique_rows[rev(order(  unique_rows[, ncol(unique_rows)-1])), ]
   # Result
   pec=rep(0, length(unique_rows[,1]))
   for(i in 1:length(summ))
   {
     pec=(unique_rows[,1]==summ[i])+pec
     
   }
   ms=min(length(cbind(unique_rows,as.logical(pec))[ ,1]),lim)
   
  hybrid_dist_boot_boot <- cbind(unique_rows,as.logical(pec))[1:ms,]
   }
   
   
   hybrid_dist_boot_boot  <- function(bootlvl1,lvl1,lim ) {
   
   pic=list()
   for(i in 1:length(bootlvl1))
   {
     pic[[i]]=hybrid_dist_boot(bootlvl1[[i]],bootlvl1)
   }
   combined <- do.call(rbind,pic)
   
   # Step 2: Remove duplicate rows based on the first column
   unique_rows <- combined[!duplicated(combined[, 1]), ]
   
   summ<- cycle_circord(lvl1)
   summ<-cycle_circord_label(summ)
   unique_rows <-   unique_rows[rev(order(  unique_rows[, ncol(unique_rows)-1])), ]
   # Result
   pec=rep(0, length(unique_rows[,1]))
   for(i in 1:length(summ))
   {
     pec=(unique_rows[,1]==summ[i])+pec
     
   }
   ms=min(length(cbind(unique_rows,as.logical(pec))[ ,1]),lim)
   
  hybrid_dist_boot_boot <- cbind(unique_rows,as.logical(pec))[1:ms,]
 }
 #############################
 #############################
 #####  Cycle  Distance  #####
 #############################
 #############################
 #############################
 #############################
 vectorize_circord <- function(net1) {
   c <- cycle_circord(net1)
   tax <- sort(net1$tip.label)
   
   lapply(seq_along(c), function(i) {
     cycle_i <- c[[i]]
     
     # Create position vector with NAs
     p <- rep(NA_integer_, length(tax))
     
     # Fill with group index using match (faster than %in%+which)
     for (j in seq_along(cycle_i)) {
       idx <- match(unlist(cycle_i[[j]]), tax)
       p[idx] <- j
     }
     
     list(seq_along(cycle_i), p)
   })
 }
 
 
 cycle_dist_boot <- function(lvl1, bootlvl1) {
   
   if(is.null(lvl1$reticulation)){cycle_dist_boot <- 'this is a tree'}else{
   
   f1 <- vectorize_circord(lvl1)
   
   if(length(f1)>0){
   circord <- cycle_circord(lvl1)
   data<-circord
   circord <- data[!sapply(data, function(inner_list) any(sapply(inner_list, length) == 0))]
   
   inull <-  which(sapply(circord, is.null)==T)
   circord <- circord[!sapply(circord, is.null)]
   labelorder <- cycle_circord_label(circord)
   if(length(inull)>0)
   {f1 <- f1[-inull]}
   
    
   n  <- length(lvl1$tip.label)
   
   cycscores <- matrix(0, nrow = length(bootlvl1), ncol = length(f1))
   
   # for (i in seq_along(bootlvl1)) 
     # f2 <- vectorize_circord(bootlvl1[[i]])
   #   cycscores[i, ] <- vapply(f1, function(cycle1) {
   #     # Compute distances to all cycles in f2
   #     dists <- vapply(f2, function(cycle2) {
   #       ordereddist(cycle1[[1]], cycle1[[2]], cycle2[[1]], cycle2[[2]])
   #     }, numeric(1))
   #     min(dists)
   #   }, numeric(1))
   # }
   # 
   # cycscores=matrix(rep(0, length(f1)*length(bootlvl1)), ncol =length(f1) )
   
   for(i in 1:length(bootlvl1))
   {
     h <-  bootlvl1[[i]]$reticulation
     
     
     net2 =  bootlvl1[[i]] 
     
     
     for(j in 1:length(f1))
     {
       zip=c()
       if(!is.null(dim(h)[1])){
         f2=vectorize_circord(net2)
       for(k in 1:length(f2))
       {
         zip[k] =  ordereddist(f1[[j]][[1]],f1[[j]][[2]],f2[[k]][[1]],f2[[k]][[2]])        
       }
       
       cycscores[i,j] = min(zip)
       }else{cycscores[i,j] =1  }
     }
     
   }
  cycledist <- list( round(1 - (colMeans(cycscores) / n),3),round( apply(cycscores, 2, function(x) {return(sum(x==0)) })/length(bootlvl1),3) )  
  # Convert to matrix (each list element becomes a row)
  cycledist <- cbind(labelorder,t(do.call(rbind, cycledist)))
  if(is.vector(cycledist)){ cycledist<- matrix(cycledist, ncol = 3)}
  cycledist <- `colnames<-`(cycledist,c("Circular order","TS","BP")) 
  
  count <- str_count( cycledist[,1],' / ')
  cycledist <-  cycledist[count>2,]
  if(is.vector(cycledist)){ cycledist<- matrix(cycledist, ncol = 3);cycledist <- `colnames<-`(cycledist,c("Circular order","TS","BP")) }
   return(cycledist)}else{return('No blobs')}
   }
 }
 
 #to deal with 3-cycles that may cause problems
 cycle_circordv2 <- function(lvl1) {
   taxa <- lvl1$tip.label
   il1  <- igraph::as.igraph(lvl1)
   
   # non-bridge edges
   cut_edges     <- igraph::bridges(il1)
   non_cut_edges <- base::setdiff(E(il1), cut_edges)
   
   # nodes with indegree = 2
   indeg             <- igraph::degree(il1, mode = "in")
   nodes_indegree_2  <- V(il1)[indeg == 2]
   
   # subgraph of non-bridge edges
   g_sub <- igraph::as_undirected(
     igraph::subgraph_from_edges(il1, eids = non_cut_edges, delete.vertices = TRUE),
     mode = "collapse"
   )
   comp <- igraph::components(g_sub)$membership
   
   # process each indegree-2 node
   data = lapply(seq_along(nodes_indegree_2), function(i) {
     cycs  <- names(comp[comp == i])
     cycle <- igraph::induced_subgraph(g_sub, cycs)
     
     girth_info <- igraph::girth(cycle, circle = TRUE)
     if (girth_info$girth == Inf) return(NULL)
     
     lapply(seq_along(girth_info$circle), function(j) {
       g3     <- igraph::delete_vertices(il1, names(girth_info$circle[j]))
       comp3  <- igraph::components(g3)$membership
       index  <- base::setdiff(c(1, 2), unique(comp3[names(girth_info$circle)[-j]]))
       intersect(names(comp3[comp3 == index]), taxa)
     })
   })
   
   if(length((which(sapply(data, function(inner_list) any(sapply(inner_list, length) == 0)))))>0){
     for( i in 1:(which(sapply(data, function(inner_list) any(sapply(inner_list, length) == 0)))) )
     {
       
       data[[i]] <-   data[[i]][lengths(data[[i]]) > 0]
     }
   }
   circord <- data 
   
   
   return(circord)
 }
 
 cycle_circord <- function(lvl1,cap=4) {
   taxa <- lvl1$tip.label
   il1  <- igraph::as.igraph(lvl1)
   
 
   cut_edges     <- igraph::bridges(il1)
   non_cut_edges <- base::setdiff(E(il1), cut_edges)
   
 
   indeg             <- igraph::degree(il1, mode = "in")
   nodes_indegree_2  <- V(il1)[indeg == 2]
   
   
   g_sub <- igraph::as_undirected(
     igraph::subgraph_from_edges(il1, eids = non_cut_edges, delete.vertices = TRUE),
     mode = "collapse"
   )
   comp <- igraph::components(g_sub)$membership
   
   
   data = lapply(seq_along(nodes_indegree_2), function(i) {
     cycs  <- names(comp[comp == i])
     cycle <- igraph::induced_subgraph(g_sub, cycs)
     
     girth_info <- igraph::girth(cycle, circle = TRUE)
     if (girth_info$girth == Inf) return(NULL)
     
     lapply(seq_along(girth_info$circle), function(j) {
       g3     <- igraph::delete_vertices(il1, names(girth_info$circle[j]))
       comp3  <- igraph::components(g3)$membership
       index  <- base::setdiff(c(1, 2), unique(comp3[names(girth_info$circle)[-j]]))
       intersect(names(comp3[comp3 == index]), taxa)
     })
   })
   
 if(length((which(sapply(data, function(inner_list) any(sapply(inner_list, length) == 0)))))>0){
    for( i in 1:(which(sapply(data, function(inner_list) any(sapply(inner_list, length) == 0)))) )
    {
 
      data[[i]] <-   data[[i]][lengths(data[[i]]) > 0]
    }
 }
     circord <- data[!sapply(data, function(inner_list)  length(inner_list ) <cap)]
   
    
   return(circord)
 }
 
 
 ordereddist <- function(o1, p1, o2, p2) {
 
   w1= c(get_cyclic_permutations(o1),get_cyclic_permutations(rev(o1)))
   w2= c(get_cyclic_permutations(o2),get_cyclic_permutations(rev(o2)))
   ss=1
   res=c()
   for(ii in 1:length(w1))
   {
     for(jj in 1:length(w2))
     { 
       tab <<- matrix(NA, nrow = length(o1), ncol = length(o2))
       res[ss] <- ordaux(1, w1[[ii]], p1, 1, w2[[jj]], p2)
       ss= ss+1
     }
   }
   return(( length(p1)  - max(res )))
 }
 
 # Recursive helper with memoization
 ordaux <- function(i, o1, p1, j, o2, p2) {
   if (i > length(o1) || j > length(o2)) {
     return(0)
   }
   
   if (!is.na(tab[i, j])) {
     return(tab[i, j])
   }
   
   s11 <- shared(o1[i], p1, o2[j], p2) + ordaux(i + 1, o1, p1, j + 1, o2, p2)
   s10 <- ordaux(i + 1, o1, p1, j, o2, p2)
   s01 <- ordaux(i, o1, p1, j + 1, o2, p2)
   
   max_val <- max(s11, s10, s01)
   tab[i, j] <<- max_val
   return(max_val)
 }
 
 # Count taxa shared by c1 in p1 and c2 in p2
 shared <- function(c1, p1, c2, p2) {
   sum(p1 == c1 & p2 == c2)
 }
 
 get_cyclic_permutations <- function(vec) {
   n <- length(vec)
   lapply(0:(n - 1), function(i) c(tail(vec, i), head(vec, n - i)))
 }
 
 #############################
 #############################
 #####  edge   Distance  #####
 #############################
 #############################
 
 split_matrix <- function(tree) {
   if (!requireNamespace("ape", quietly = TRUE)) {
     stop("Package 'ape' is required")
   }
   
   tree <- ape::unroot(tree) 
   taxa <- sort(tree$tip.label)   
   tree <- root(tree,outgroup =  taxa[1] )# Treat tree as unrooted
      # Columns in alphabetical order
   splits <- prop.part(tree)
     
   # Initialize matrix
   mat <- matrix(
     0,
     nrow = length(splits),
     ncol = length(taxa),
     dimnames = list(
       paste0("split_", seq_along(splits)),
       taxa
     )
   )
   
   for (i in seq_along(splits)) {
     split_taxa <- tree$tip.label[splits[[i]]]
     
     # FIX: always use the smaller side of the split
     if (length(split_taxa) > length(taxa) / 2) {
       split_taxa <- setdiff(taxa, split_taxa)
     }
     
     # Fill matrix
     mat[i, taxa %in% split_taxa] <- 1
   }
   
   return(mat)
 }
 
 
 split_distance <- function(s, t) {
   if (length(s) != length(t))
     stop("Splits must have same length")
   
   if (!all(s %in% c(0,1)) || !all(t %in% c(0,1)))
     stop("Splits must be binary (0/1)")
   
   d1 <- sum(abs(s - t))        # direct
   d2 <- sum(abs(s - (1 - t)))  # flipped
   
   min(d1, d2)
 }
 
 tree_split_distances <- function(tree,M1, tree1) {
   
   tree=unroot.phylo(tree)
   
   tree   <- drop.tip(tree, "[ROOT]")
   tree1=unroot.phylo(tree1)
   taxa <- sort(tree$tip.label)  
   
   
   M2 <- split_matrix(tree1)[-1,]
   
   if (!identical(colnames(M1), colnames(M2))) {
     stop("Trees must have identical taxa sets")
   }
   
   distances <-  1-  apply(M1, 1, function(s) {
   min(apply(M2, 1, function(t) split_distance(s, t)))
   })/floor(length(taxa)/2)  
   
   tree_split_distances <-data.frame(
     split = rownames(M1),
     distance = distances,
     row.names = NULL
   )
 }
 
edge_dist_boot  <- function(tree, boot_trees) {
  M1 <- split_matrix(tree)[-1,]
  TS =rep(0,(dim(split_matrix(tree)[-1,])[1]))
  BP =rep(0,(dim(split_matrix(tree)[-1,])[1]))
  for(i in 1:length(boot_trees))
  {
   TS =   TS+tree_split_distances(tree,M1,boot_trees[[i]])[,2] 
    sd= tree_split_distances(tree,M1,boot_trees[[i]])[,2]==1
    
    BP =  BP+ sd
  }
   
 
  names  <- sapply(1:nrow(M1), function(i) {
    row <- M1[i, ]
    target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })
   support =  rbind( TS,BP) 
    colnames(support)<-names
  edge_dist_boot  <-   support  
} 
merge_outputs <- function(out1, out2, out3) {
  out1 <- as.data.frame(out1)
  out2 <- as.data.frame(out2)
  out3 <- as.data.frame(out3)
  
  all_cols <- unique(c(colnames(out1), colnames(out2), colnames(out3)))
  
  align <- function(out) {
    missing <- all_cols[!all_cols %in% colnames(out)]
    if (length(missing) > 0) {
      out[, missing] <- NA
    }
    out[, all_cols]
  }
  
  result <- rbind(align(out1), align(out2), align(out3))
  rownames(result) <- c(paste0(rownames(out1), "_out1"),
                        paste0(rownames(out2), "_out2"),
                        paste0(rownames(out3), "_out3"))
  result
}
 #############################
 #############################
 #############################
 #############################
 #####  Blob   Distance  #####
 #############################
 #############################
 

 obtain_tree_partition = function(Tr)
 {
   names=Tr$tip.label
   
   iTr=as.igraph(Tr, directed = FALSE)
   graph.attributes(iTr)
   highdeg = V(iTr)[which(igraph::degree(iTr)>2)]
   nodes = list()
   if(sum(igraph::degree(iTr)>2)>0){
     k=1
     for(v in highdeg)
     {
       df= as.data.frame(components(delete_vertices(iTr, v))$membership)
       df <- tibble::rownames_to_column(df, "VALUE")
       nodes[[k]] =  df[df$VALUE%in%names,]
       k=k+1
     }
   }
   
   
   
   PRT=list()
   for(i in 1:length(nodes))
   {
     l= nodes[[i]]
     PRT[[i]] = list()
     for(j in 1:max(l[,2]))
     {
       PRT[[i]][[j]] = sort( l[which(l[,2]==j),1] )
     }
   }
   
   PRT2=list()
   nodes2=list()
   k=1
   for(i in 1:length(nodes))
   {
     if(length(unlist(unique(nodes[[i]][2])))>3)
     {
       PRT2[[k]]= PRT[[i]]
       nodes2[[k]] = nodes[[i]]
       k=k+1
     }
     
   }
   
   obtain_tree_partition = list(PRT2,nodes2)
 }
 
 TF_dist <- function(x, y) {
   xn <- max(x[,2])
   yn <- max(y[,2])
   
   # Precompute row sets for x and y
   x_sets <- lapply(1:xn, function(i) x[x[,2] == i, 1])
   y_sets <- lapply(1:yn, function(j) y[y[,2] == j, 1])
   
   # Initialize cost matrix
   C <- matrix(0L, nrow = xn, ncol = yn)
   
   for (i in seq_len(xn)) {
     for (j in seq_len(yn)) {
       inter_len <- length(intersect(x_sets[[i]], y_sets[[j]]))
       C[i, j] <- inter_len
     }
   }
   
   # Hungarian algorithm on -C (maximize overlap)
   cost <- HungarianSolver(-C)$cost
   
   # Total TF distance
   nrow(x) + cost
 }
 
 DeltaB<-function(x)
 {
   n <- dim(x)[1]
   k <- max(x[,2])
   
   mm <- min(0, min(sort(table(x[,2]))-(1:k)) )
   DeltaB <- n -k -mm
 }
  
 TI <- function(x,Tree)
 {
   
   treepart <- obtain_tree_partition(Tree)[[2]]
   disc<-c()
   n <- dim(x)[1]
   if(length(treepart)>0)
   {
     for(i in 1:length(treepart))
     {
       disc[i] <- TF_dist(x,treepart[[i]])
     }
     
     TI <- min(disc)
   }else{TI<-n}
   
 }
   
 TS <- function(x, Tree) {
   n <- nrow(x)
   1 - TI(x, Tree) / n
 }
 
 
 blob_dist_boot  <- function(Tree, boot_trees) {
   prT <- obtain_tree_partition(Tree)[[2]]
   

   if(length(prT)>0)
   {
   prT <-  lapply(prT, as.data.frame) 
   suppts   <- numeric(length(prT))   # average TS support
   supptb <- numeric(length(prT))   # count of exact TS == 1
   
   for (i in seq_along(prT)) {
     for (k in seq_along(boot_trees)) {
       ts_val <- TS(prT[[i]], boot_trees[[k]])
       suppts[i]   <- suppts[i]   + ts_val
       supptb[i] <- supptb[i] + as.integer(ts_val == 1)
     }
   }
   
   suppts   <- suppts   / length(boot_trees)
   supptb <- supptb / length(boot_trees)   # normalize to proportion if you want
   
   supp = round(cbind(suppts,supptb),3)
 
   blobtable <- summary_ToB(Tree, t(supp)) 
    
    count <- str_count( blobtable[,1],'/')
    blobtable <-  blobtable[count>2,]
    if(is.vector(blobtable)){blobtable = matrix(blobtable, ncol = 3)}
    blobtable<-  `colnames<-`(blobtable,c('Blob',"TS",'BP'))
  
 return(blobtable)
   }else{'No blobs'}
 }
 
  
 summary_ToB <- function(Tree, supp = FALSE) {
   prT <- obtain_tree_partition(Tree)[[2]]
   nodes=c()
   for(i in 1:length(prT))
   {
     x<-as.matrix(prT[[i]])
     xx=as.numeric(x[,2])
     part<-c()
     k <- max(xx)
     for(j in 1:k)
     {
       part<-c(part,  paste0(sort(x[,1][which(xx==j)]) ,collapse = '*'))
     }
     
     nodes[i]<-paste0(sort(part),collapse = ' / ')
   }
   if(   all(is.logical(supp)) )
   {
     summary_toB<- t(rbind(nodes))
   }else {summary_toB<- t(rbind(nodes,supp))
   if(is.vector(summary_toB)){summary_toB = matrix(summary_toB, ncol = 3)}
   summary_toB <-  `colnames<-`(summary_toB,c('Blob',"TS",'BP'))
   }
   summary_toB
 }
 
 
 blob_dist_boot_boot  <- function(boot_trees, Tree,lim) {
  
   pic=list()
   for(i in 1:length(boot_trees))
   {
     pic[[i]]=blob_dist_boot(boot_trees[[i]],boot_trees)
   }
   combined <- do.call(rbind,pic)
   
   # Step 2: Remove duplicate rows based on the first column
   unique_rows <- combined[!duplicated(combined[, 1]), ]
   
   summ<-summary_ToB(Tree, supp = FALSE)
   unique_rows <-   unique_rows[rev(order(  unique_rows[, ncol(unique_rows)-1])), ]
   # Result
   pec=rep(0, length(unique_rows[,1]))
   for(i in 1:length(summ))
   {
     pec=(unique_rows[,1]==summ[i])+pec
     
   }
   ms=min(length(cbind(unique_rows,as.logical(pec))[ ,1]),lim)
   
   blob_dist_boot_boot <- cbind(unique_rows,as.logical(pec))[1:ms,]
 }
 #############################
 #############################
 ##### Auxiliar functions ####
 #############################
 #############################
 #############################
 #############################
 read.evonets <- function(file)
 {
   L <- readLines(file)
   L <- L[!L=="fullResMax exceeded; full resolutions not produced"]
   nets <- list()
   for(i in 1:length(L))
   {
     
     nets[[i]]=read.evonet(text = L[i])
     if(dim(nets[[i]]$reticulation)[1]==0)
     {
       nets[[i]]=read.tree(text = L[i])
     }
   }
   return(nets)
 }
 
 read.nets <- function(file)
 {
   L <- readLines(file)
   L <- L[!L=="fullResMax exceeded; full resolutions not produced"]
   nets <- list()
   for(i in 1:length(L))
   {
     
     nets[[i]]=read.net(text = L[i])
     if(dim(nets[[i]]$reticulation)[1]==0)
     {
       nets[[i]]=read.tree(text = L[i])
     }
   }
   return(nets)
 }
 
 Trees_of_Blobs <- function(networks)
 {
 
   return(lapply(networks, treeOfBlobs))
 }
 
 Boostrap_Tinnik <- function(genedata, sample_size = F, num_trees = 100, omit = FALSE, epsilon = 0, test = "T3", alpha = 0.05, beta = 0.95, treemethod = fastme.bal, delta = 2, taxanames = NULL)
{
   if (is.null(sample_size) || isFALSE(sample_size)) {
     sample_size <- length(genedata)
   }
   
   # Generate bootstrap trees
   boot_trees <- lapply(seq_len(num_trees), function(i) {
     gts <- sample(genedata, sample_size, replace = TRUE)
     TINNIK(gts, omit, epsilon, test, alpha, beta, treemethod, delta, taxanames, plot = FALSE)$ToB
   })
   
   boot_trees
 }
 
 plot_evonet_with_gamma <- function(enewick,
                                    digits = 2,
                                    tip_text_size = 3,
                                    gamma_text_size = 3,
                                    tick_text_size = 3,
                                    plot_title = NULL,
                                    scale_pad = 0.06,
                                    scale_gap = 0.12,
                                    tick_len = 0.03,
                                    pad_x_frac = 0.25,
                                    pad_y_frac = 0.15,
                                    gamma_nudge = 0.05,
                                    show_scale = TRUE) {  # <-- new parameter
   stopifnot(requireNamespace("ape", quietly = TRUE))
   stopifnot(requireNamespace("ggplot2", quietly = TRUE))
   if (!requireNamespace("ggnetworx", quietly = TRUE) &&
       !requireNamespace("tanggle", quietly = TRUE)) {
     stop("Need 'ggnetworx' (GitHub) or 'tanggle' (Bioconductor) for ggevonet().")
   }
   ggevonet_fun <- if (requireNamespace("ggnetworx", quietly = TRUE)) ggnetworx::ggevonet else tanggle::ggevonet
   
   # 1) Base plot
   net <- ape::read.evonet(text = enewick)
   p <- ggevonet_fun(net, ladderize = FALSE) +
     ggtree::geom_tiplab(size = tip_text_size) +
     ggplot2::theme(legend.position = "none")
   
   df <- p$data
   
   # --- compute depth and diameter ---
   xmin <- min(df$x, na.rm = TRUE)
   xmax <- max(df$x, na.rm = TRUE)
   total_depth <- xmax - xmin
   diameter <- 2 * total_depth
   
   if (!is.null(plot_title)) {
     # full_title <- paste0(plot_title,
                          # "  |  Network diameter in cu: ",
                          # sprintf(paste0("%.", digits, "f"), diameter))
     
     full_title <- paste0(plot_title )
     p <- p + ggplot2::ggtitle(full_title) +
       ggplot2::theme(plot.title.position = "plot")
   }
   
   # --- parse gamma values ---
   rx <- gregexpr("#H\\d+:[^():;,]*::([0-9eE+\\-\\.]+)", enewick, perl = TRUE)
   hits <- regmatches(enewick, rx)[[1]]
   if (length(hits) > 0) {
     get_hid   <- function(s) sub(".*?(#H\\d+):.*", "\\1", s, perl = TRUE)
     get_gamma <- function(s) as.numeric(sub(".*?::([0-9eE+\\-\\.]+).*", "\\1", s, perl = TRUE))
     gtab <- data.frame(
       hid = vapply(hits, get_hid, character(1)),
       gamma = vapply(hits, get_gamma, numeric(1)),
       stringsAsFactors = FALSE
     )
     gtab <- gtab[!duplicated(gtab$hid), , drop = FALSE]
   } else {
     gtab <- data.frame(hid = character(), gamma = numeric())
   }
   
   # --- gamma labels with offset ---
   parents <- df[, c("node","x","y")]
   names(parents) <- c("parent","parent.x","parent.y")
   df2 <- merge(df, parents, by = "parent", all.x = TRUE, sort = FALSE)
   
   extract_h <- function(lbl) {
     if (is.null(lbl)) return(NA_character_)
     m <- regexpr("#H\\d+", lbl, perl = TRUE)
     ifelse(m > 0, regmatches(lbl, m), NA_character_)
   }
   df2$hid <- vapply(df2$label, extract_h, character(1))
   ann <- merge(
     df2[!is.na(df2$hid), c("node","parent","x","y","parent.x","parent.y","hid")],
     transform(gtab, lab = sprintf(paste0("%.", digits, "f"), gamma)),
     by = "hid", all.x = TRUE, all.y = FALSE
   )
   ann <- ann[!is.na(ann$lab), , drop = FALSE]
   
   if (nrow(ann) > 0) {
     ann$dx <- ann$x - ann$parent.x
     ann$dy <- ann$y - ann$parent.y
     len <- sqrt(ann$dx^2 + ann$dy^2)
     ann$offx <- ifelse(len > 0, -ann$dy/len, 0)
     ann$offy <- ifelse(len > 0,  ann$dx/len, 0)
     ann$mx <- (ann$x + ann$parent.x)/2 + gamma_nudge * ann$offx
     ann$my <- (ann$y + ann$parent.y)/2 + gamma_nudge * ann$offy
   }
   
   # --- coord limits (adjusted depending on show_scale) ---
   ymin <- min(df$y, na.rm = TRUE)
   ymax <- max(df$y, na.rm = TRUE)
   xpad <- pad_x_frac * (xmax - xmin)
   ypad <- pad_y_frac * max(1e-8, (ymax - ymin))
   
   if (show_scale) {
     D_seq  <- seq(0, floor(total_depth), by = 1)
     ybase  <- ymin - (scale_pad + scale_gap)
     tick_x <- xmax - D_seq
     scale_df <- data.frame(x = c(xmax, xmin), y = c(ybase, ybase))
     ticks_df <- data.frame(
       x = tick_x, xend = tick_x,
       y = ybase - tick_len/2, yend = ybase + tick_len/2,
       lab = D_seq
     )
     
     p <- p +
       ggplot2::coord_cartesian(
         xlim = c(xmin - xpad, xmax + xpad),
         ylim = c(ybase - tick_len - scale_pad/2 - ypad, ymax + ypad),
         expand = FALSE
       ) +
       ggplot2::geom_segment(
         data = scale_df[1, , drop = FALSE],
         ggplot2::aes(x = x, xend = scale_df[2,1], y = y, yend = y),
         linewidth = 0.4, inherit.aes = FALSE
       ) +
       ggplot2::geom_segment(
         data = ticks_df,
         ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
         linewidth = 0.3, inherit.aes = FALSE
       ) +
       ggplot2::geom_text(
         data = ticks_df,
         ggplot2::aes(x = x, y = yend + tick_len*0.8, label = lab),
         size = tick_text_size, vjust = 0, inherit.aes = FALSE
       )
   } else {
     p <- p +
       ggplot2::coord_cartesian(
         xlim = c(xmin - xpad, xmax + xpad),
         ylim = c(ymin - ypad, ymax + ypad),
         expand = FALSE
       )
   }
   
   # gamma labels
   if (nrow(ann) > 0) {
     p <- p + ggplot2::geom_text(
       data = ann,
       ggplot2::aes(x = mx, y = my, label = lab),
       color = "red", size = gamma_text_size, inherit.aes = FALSE
     )
   }
   
   p
 }
 
 NewickFromBipart <- function(filename) {
   lines <- readLines(filename)
   
   bipart_idx <- which(lines == "BIPARTITIONS non redundant with blobs")
   taxa_idx <- which(lines == "ENCODING_OF_TAXA")
   blob_idx <- which(lines == "BLOBS given as cycles, where 0 is the hybrid clade")
   
   taxa_section <- lines[(taxa_idx + 2):(blob_idx - 2)]
   taxa_names <- sub(".*\\s([^)]*)\\)", "\\1", taxa_section)
   
   bipartitions <- lines[(bipart_idx + 2):(length(lines) - 1)]
   blobs <- lines[(blob_idx + 2):(bipart_idx - 2)]
   
   num_taxa <- length(taxa_names)
   
   tree <- rmtree(2, n = num_taxa, tip.label = 1:num_taxa)
   q_table <- quartetTable(tree)
   
   q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] <- 
     q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] * 0 + 50
   
   selected_rows <- c()
   attributes(q_table)$dimnames[[2]] <- c(taxa_names, "12|34", "13|24", "14|23", "1234")
   
   for (i in 1:length(bipartitions)) {
     seq_str <- gsub(" ", "", bipartitions[i])
     bits <- strsplit(seq_str, "")[[1]]
     
     pos_ones <- which(bits == "1")
     pos_zeroes <- which(bits == "0")
     
     comb_ones <- combn(pos_ones, 2)
     comb_zeroes <- combn(pos_zeroes, 2)
     
     for (ii in 1:(dim(comb_ones)[2])) {
       xx <- sort(comb_ones[, ii])
       
       for (jj in 1:(dim(comb_zeroes)[2])) {
         yy <- sort(comb_zeroes[, jj])
         
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[c(xx, yy)] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (xx[2] < yy[1] || yy[2] < xx[1]) {
           tag <- c(100, 25, 25, 0)
         } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                    (yy[1] < xx[1] && xx[2] < yy[2])) {
           tag <- c(25, 25, 100, 0)
         } else {
           tag <- c(25, 100, 25, 0)
         }
         
         q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
       }
     }
   }
   
   for (i in 1:length(blobs)) {
     seq_str <- gsub(" ", "", blobs[i])
     bits <- as.numeric(strsplit(seq_str, "")[[1]])
     csize <- max(bits) + 1
     
     pos_zeroes <- which(bits == 0)
     pos_list <- list(pos_zeroes)
     
     for (z in 2:csize) {
       pos_list[[z]] <- which(as.numeric(bits) == z - 1)
     }
     
     combs <- combn(1:csize, 4)
     
     for (z in 1:dim(combs)[2]) {
       desk <- sort(combs[, z])
       grid <- expand.grid(
         pos_list[desk[1]][[1]],
         pos_list[desk[2]][[1]],
         pos_list[desk[3]][[1]],
         pos_list[desk[4]][[1]]
       )
       
       for (zz in 1:(dim(grid)[1])) {
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[unlist(grid[zz,])] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (1 %in% desk) {
           tag <- c(50, 50, 50, 0)
           
           for (zzz in 1:2) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(25, 0, 0, 0) + tag
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) ||
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(0, 0, 25, 0) + tag
             } else {
               tag <- c(0, 25, 0, 0) + tag
             }
           }
           
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
           
         } else {
           for (zzz in 1) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(100, 25, 25, 0)
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(25, 25, 100, 0)
             } else {
               tag <- c(25, 100, 25, 0)
             }
           }
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
         }
       }
     }
     
     for (sa in 1:length(pos_list)) {
       if (length(pos_list[[sa]]) > 1) {
         dd <- combn(pos_list[[sa]], 2)
         
         combined <- unlist(pos_list)
         combined <- base::setdiff(combined, pos_list[[sa]])
         dds <- combn(combined, 2)
         
         for (isd in 1:(dim(dd)[2])) {
           xx <- sort(as.vector(dd[, isd]))
           
           for (isds in 1:(dim(dds)[2])) {
             yy <- sort(as.vector(dds[, isds]))
             
             target <- c(rep(0, num_taxa), 50, 50, 50, 0)
             target[c(xx, yy)] <- 1
             
             drow <- which(apply(q_table, 1, function(row) all(row == target)))
             
             if (length(drow) > 0) {
               if (xx[2] < yy[1] || yy[2] < xx[1]) {
                 tag <- c(100, 25, 25, 0)
               } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                          (yy[1] < xx[1] && xx[2] < yy[2])) {
                 tag <- c(25, 25, 100, 0)
               } else {
                 tag <- c(25, 100, 25, 0)
               }
               
               q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
             }
           }
         }
       }
     }
   }
   
   tob <- TINNIK(q_table)
   ss <- resolveLevel1(tob$ToB, tob$pTable, alpha = 0.05, beta = 0.95)
   
   return(ss[[1]][[1]])
 }
 
 
 
 NewickFromBipart_parts <- function(taxa_names,bipartitions,blobs) {
 
   num_taxa <- length(taxa_names)
   
   tree <- rmtree(2, n = num_taxa, tip.label = 1:num_taxa)
   q_table <- quartetTable(tree)
   
   q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] <- 
     q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] * 0 + 50
   
   selected_rows <- c()
   attributes(q_table)$dimnames[[2]] <- c(taxa_names, "12|34", "13|24", "14|23", "1234")
   
   for (i in 1:length(bipartitions)) {
     seq_str <- gsub(" ", "", bipartitions[i])
     bits <- strsplit(seq_str, "")[[1]]
     
     pos_ones <- which(bits == "1")
     pos_zeroes <- which(bits == "0")
     
     comb_ones <- combn(pos_ones, 2)
     comb_zeroes <- combn(pos_zeroes, 2)
     
     for (ii in 1:(dim(comb_ones)[2])) {
       xx <- sort(comb_ones[, ii])
       
       for (jj in 1:(dim(comb_zeroes)[2])) {
         yy <- sort(comb_zeroes[, jj])
         
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[c(xx, yy)] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (xx[2] < yy[1] || yy[2] < xx[1]) {
           tag <- c(100, 25, 25, 0)
         } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                    (yy[1] < xx[1] && xx[2] < yy[2])) {
           tag <- c(25, 25, 100, 0)
         } else {
           tag <- c(25, 100, 25, 0)
         }
         
         q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
       }
     }
   }
   
   for (i in 1:length(blobs)) {
     seq_str <- gsub(" ", "", blobs[i])
     bits <- as.numeric(strsplit(seq_str, "")[[1]])
     csize <- max(bits) + 1
     
     pos_zeroes <- which(bits == 0)
     pos_list <- list(pos_zeroes)
     
     for (z in 2:csize) {
       pos_list[[z]] <- which(as.numeric(bits) == z - 1)
     }
     
     combs <- combn(1:csize, 4)
     
     for (z in 1:dim(combs)[2]) {
       desk <- sort(combs[, z])
       grid <- expand.grid(
         pos_list[desk[1]][[1]],
         pos_list[desk[2]][[1]],
         pos_list[desk[3]][[1]],
         pos_list[desk[4]][[1]]
       )
       
       for (zz in 1:(dim(grid)[1])) {
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[unlist(grid[zz,])] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (1 %in% desk) {
           tag <- c(50, 50, 50, 0)
           
           for (zzz in 1:2) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(25, 0, 0, 0) + tag
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) ||
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(0, 0, 25, 0) + tag
             } else {
               tag <- c(0, 25, 0, 0) + tag
             }
           }
           
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
           
         } else {
           for (zzz in 1) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(100, 25, 25, 0)
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(25, 25, 100, 0)
             } else {
               tag <- c(25, 100, 25, 0)
             }
           }
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
         }
       }
     }
     
     for (sa in 1:length(pos_list)) {
       if (length(pos_list[[sa]]) > 1) {
         dd <- combn(pos_list[[sa]], 2)
         
         combined <- unlist(pos_list)
         combined <- base::setdiff(combined, pos_list[[sa]])
         dds <- combn(combined, 2)
         
         for (isd in 1:(dim(dd)[2])) {
           xx <- sort(as.vector(dd[, isd]))
           
           for (isds in 1:(dim(dds)[2])) {
             yy <- sort(as.vector(dds[, isds]))
             
             target <- c(rep(0, num_taxa), 50, 50, 50, 0)
             target[c(xx, yy)] <- 1
             
             drow <- which(apply(q_table, 1, function(row) all(row == target)))
             
             if (length(drow) > 0) {
               if (xx[2] < yy[1] || yy[2] < xx[1]) {
                 tag <- c(100, 25, 25, 0)
               } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                          (yy[1] < xx[1] && xx[2] < yy[2])) {
                 tag <- c(25, 25, 100, 0)
               } else {
                 tag <- c(25, 100, 25, 0)
               }
               
               q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
             }
           }
         }
       }
     }
   }
   
   tob <- TINNIK(q_table)
   ss <- resolveLevel1(tob$ToB, tob$pTable, alpha = 0.05, beta = 0.95  )
   sos <- NANUQ(tob$ToB, tob$pTable, alpha = 0.05, beta = 0.95  )
   
   return(ss[[1]][[1]])
 }
 
 NewickFromBipart_partsv2 <- function(taxa_names,bipartitions,blobs) {
   
   num_taxa <- length(taxa_names)
   
   tree <- rmtree(2, n = num_taxa, tip.label = 1:num_taxa)
   q_table <- quartetTable(tree)
   
   q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] <- 
     q_table[, (dim(q_table)[2] - 3):(dim(q_table)[2] - 1)] * 0 + 50
   
   selected_rows <- c()
   attributes(q_table)$dimnames[[2]] <- c(taxa_names, "12|34", "13|24", "14|23", "1234")
   
   for (i in 1:length(bipartitions)) {
     seq_str <- gsub(" ", "", bipartitions[i])
     bits <- strsplit(seq_str, "")[[1]]
     
     pos_ones <- which(bits == "1")
     pos_zeroes <- which(bits == "0")
     
     comb_ones <- combn(pos_ones, 2)
     comb_zeroes <- combn(pos_zeroes, 2)
     
     for (ii in 1:(dim(comb_ones)[2])) {
       xx <- sort(comb_ones[, ii])
       
       for (jj in 1:(dim(comb_zeroes)[2])) {
         yy <- sort(comb_zeroes[, jj])
         
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[c(xx, yy)] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (xx[2] < yy[1] || yy[2] < xx[1]) {
           tag <- c(100, 25, 25, 0)
         } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                    (yy[1] < xx[1] && xx[2] < yy[2])) {
           tag <- c(25, 25, 100, 0)
         } else {
           tag <- c(25, 100, 25, 0)
         }
         
         q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
       }
     }
   }
   
   if(length(blobs)>0)
   {
   for (i in 1:length(blobs)) {
     seq_str <- gsub(" ", "", blobs[i])
     bits <- as.numeric(strsplit(seq_str, "")[[1]])
     csize <- max(bits) + 1
     
     pos_zeroes <- which(bits == 0)
     pos_list <- list(pos_zeroes)
     
     for (z in 2:csize) {
       pos_list[[z]] <- which(as.numeric(bits) == z - 1)
     }
     
     combs <- combn(1:csize, 4)
     
     for (z in 1:dim(combs)[2]) {
       desk <- sort(combs[, z])
       grid <- expand.grid(
         pos_list[desk[1]][[1]],
         pos_list[desk[2]][[1]],
         pos_list[desk[3]][[1]],
         pos_list[desk[4]][[1]]
       )
       
       for (zz in 1:(dim(grid)[1])) {
         target <- c(rep(0, num_taxa), 50, 50, 50, 0)
         target[unlist(grid[zz,])] <- 1
         
         drow <- which(apply(q_table, 1, function(row) all(row == target)))
         selected_rows <- c(selected_rows, drow)
         
         if (1 %in% desk) {
           tag <- c(50, 50, 50, 0)
           
           for (zzz in 1:2) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(25, 0, 0, 0) + tag
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) ||
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(0, 0, 25, 0) + tag
             } else {
               tag <- c(0, 25, 0, 0) + tag
             }
           }
           
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
           
         } else {
           for (zzz in 1) {
             xx <- sort(unlist(grid[zz,][c(1, 2 * zzz)]))
             yy <- sort(unlist(base::setdiff(grid[zz,], xx)))
             
             if (xx[2] < yy[1] || yy[2] < xx[1]) {
               tag <- c(100, 25, 25, 0)
             } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                        (yy[1] < xx[1] && xx[2] < yy[2])) {
               tag <- c(25, 25, 100, 0)
             } else {
               tag <- c(25, 100, 25, 0)
             }
           }
           q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
         }
       }
     }
     
     for (sa in 1:length(pos_list)) {
       if (length(pos_list[[sa]]) > 1) {
         dd <- combn(pos_list[[sa]], 2)
         
         combined <- unlist(pos_list)
         combined <- base::setdiff(combined, pos_list[[sa]])
         dds <- combn(combined, 2)
         
         for (isd in 1:(dim(dd)[2])) {
           xx <- sort(as.vector(dd[, isd]))
           
           for (isds in 1:(dim(dds)[2])) {
             yy <- sort(as.vector(dds[, isds]))
             
             target <- c(rep(0, num_taxa), 50, 50, 50, 0)
             target[c(xx, yy)] <- 1
             
             drow <- which(apply(q_table, 1, function(row) all(row == target)))
             
             if (length(drow) > 0) {
               if (xx[2] < yy[1] || yy[2] < xx[1]) {
                 tag <- c(100, 25, 25, 0)
               } else if ((xx[1] < yy[1] && yy[2] < xx[2]) || 
                          (yy[1] < xx[1] && xx[2] < yy[2])) {
                 tag <- c(25, 25, 100, 0)
               } else {
                 tag <- c(25, 100, 25, 0)
               }
               
               q_table[drow, (dim(q_table)[2] - 3):dim(q_table)[2]] <- tag
             }
           }
         }
       }
     }
   }
   }
   # tob <- TINNIK(q_table)
   # ss <- resolveLevel1(tob$ToB, tob$pTable, alpha = 0.05, beta = 0.95  )
   # sos <- NANUQ(tob$ToB, tob$pTable, alpha = 0.05, beta = 0.95  )
   # 
   return(q_table)
 }
 
 canonical_partition <- function(x) {
   sapply(x, function(part) {
     groups <- strsplit(part, " / ")[[1]]
     
     groups_sorted <- sapply(groups, function(g) {
       elems <- strsplit(g, "\\*")[[1]]
       paste(sort(elems), collapse="*")
     })
     
     paste(sort(groups_sorted), collapse=" / ")
   })
 }
 
 merge_three_trees <- function(mat1, mat2, mat3) {
   
   df1 <- as.data.frame(mat1, stringsAsFactors = FALSE)
   df2 <- as.data.frame(mat2, stringsAsFactors = FALSE)
   df3 <- as.data.frame(mat3, stringsAsFactors = FALSE)
   
   df1$canon <- canonical_partition(df1[[1]])
   df2$canon <- canonical_partition(df2[[1]])
   df3$canon <- canonical_partition(df3[[1]])
   
   df1_sub <- df1[, c("canon", "TS", "BP")]
   colnames(df1_sub) <- c("canon", "TS_1", "BP_1")
   
   df2_sub <- df2[, c("canon", "TS", "BP")]
   colnames(df2_sub) <- c("canon", "TS_2", "BP_2")
   
   df3_sub <- df3[, c("canon", "TS", "BP")]
   colnames(df3_sub) <- c("canon", "TS_3", "BP_3")
   
   merged <- merge(df1_sub, df2_sub, by="canon", all=TRUE)
   merged <- merge(merged, df3_sub, by="canon", all=TRUE)
   
   # number of partitions
   merged$n_partitions <- sapply(strsplit(merged$canon, " / "), length)
   
   merged <- merged[, c("canon", "n_partitions",
                        setdiff(names(merged), c("canon","n_partitions")))]
   
  
   merged <- merged[!duplicated(merged$canon), ]
   
   return(merged)
 }
 
 
 merge_six_trees <- function(mat1, mat2, mat3, mat4, mat5, mat6) {
   
   df1 <- as.data.frame(mat1, stringsAsFactors = FALSE)
   df2 <- as.data.frame(mat2, stringsAsFactors = FALSE)
   df3 <- as.data.frame(mat3, stringsAsFactors = FALSE)
   df4 <- as.data.frame(mat4, stringsAsFactors = FALSE)
   df5 <- as.data.frame(mat5, stringsAsFactors = FALSE)
   df6 <- as.data.frame(mat6, stringsAsFactors = FALSE)
   
   df1$canon <- canonical_partition(df1[[1]])
   df2$canon <- canonical_partition(df2[[1]])
   df3$canon <- canonical_partition(df3[[1]])
   df4$canon <- canonical_partition(df4[[1]])
   df5$canon <- canonical_partition(df5[[1]])
   df6$canon <- canonical_partition(df6[[1]])
   
   df1_sub <- df1[, c("canon", "TS", "BP")]
   colnames(df1_sub) <- c("canon", "TS_1", "BP_1")
   
   df2_sub <- df2[, c("canon", "TS")]
   colnames(df2_sub) <- c("canon", "TS_2")
   
   df3_sub <- df3[, c("canon", "TS", "BP")]
   colnames(df3_sub) <- c("canon", "TS_3", "BP_3")
   
   df4_sub <- df4[, c("canon", "TS")]
   colnames(df4_sub) <- c("canon", "TS_4")
   
   df5_sub <- df5[, c("canon", "TS", "BP")]
   colnames(df5_sub) <- c("canon", "TS_5", "BP_5")
   
   df6_sub <- df6[, c("canon", "TS")]
   colnames(df6_sub) <- c("canon", "TS_6")
   
   merged <- merge(df1_sub, df2_sub, by="canon", all=TRUE)
   merged <- merge(merged, df3_sub, by="canon", all=TRUE)
   merged <- merge(merged, df4_sub, by="canon", all=TRUE)
   merged <- merge(merged, df5_sub, by="canon", all=TRUE)
   merged <- merge(merged, df6_sub, by="canon", all=TRUE)
   
   merged$n_partitions <- sapply(strsplit(merged$canon, " / "), length)
   
   merged <- merged[, c("canon", "n_partitions",
                        setdiff(names(merged), c("canon","n_partitions")))]
   
    
   merged <- merged[!duplicated(merged$canon), ]
   
   return(merged)
 }
 
 
 
 get_taxa_below_2cycles <- function(evonet) {
   library(igraph)
   
   # Convert evonet to igraph
   ig <- as.igraph.evonet(evonet)
   
   # Find parallel edges (2-cycles) in the igraph
   # which_multiple() identifies edges that have duplicates
   parallel_edges <- which_multiple(ig)
   
   if(sum(parallel_edges) == 0) {
     return(list(taxa = character(0), count = 0))
   }
   
   # Get edge list
   el <- as_edgelist(ig)
   
   # Find unique pairs with parallel edges
   parallel_pairs <- el[parallel_edges, , drop = FALSE]
   unique_pairs <- unique(parallel_pairs)
   
   all_taxa_below <- c()
   
   # For each parallel edge pair
   for(i in 1:nrow(unique_pairs)) {
     from_node <- unique_pairs[i, 1]
     to_node <- unique_pairs[i, 2]
     
     # Find all descendants of to_node
     descendants <- subcomponent(ig, to_node, mode = "out")
     
     # Filter for degree 1 nodes (leaves/taxa)
     degrees <- igraph::degree(ig, v = descendants, mode = "out")
     taxa_nodes <- descendants[degrees == 0]  # degree 0 out = leaves
     
     # Get labels
     taxa_labels <- V(ig)$name[taxa_nodes]
     all_taxa_below <- c(all_taxa_below, taxa_labels)
   }
   
   return(list(
     taxa = unique(all_taxa_below),
     count = length(unique(all_taxa_below))
   ))
 }
 
 # Use it
 
 ######## ######## ######## ##### 
######## Example of runs ######## 
 ######## ######## ######## ##### 
# lvl1=read.net(file ='sampled_intlvl1_3ktest6.tree')
# lvl1= bootlvl1[[100]]
# tob=treeOfBlobs(lvl1)
# bootlvl1= read.nets('sampled_bootintlvl1_3ktest6.tree'  )
# boottobs=Trees_of_Blobs(bootlvl1)
# 
# 
# Blobsupport = blob_dist_boot(tob,boottobs)
# Cyclesupport =  cycle_dist_boot( lvl1, bootlvl1)
# Hybridsupprt =  hybrid_dist_boot( lvl1, bootlvl1)
# 
# Blobsupport
# Cyclesupport
# Hybridsupprt

#  

 # Newickfrombipart<- function(filename)
 # {
 #   L = readLines(filename)
 #   bipart=which(L=="BIPARTITIONS non redundant with blobs")
 #   taxa=which(L=="ENCODING_OF_TAXA")
 #   blob=which(L=="BLOBS given as cycles, where 0 is the hybrid clade")
 #   
 #   tax = L[(taxa+2):(blob-2)]
 #   taxon <- sub(".*\\s([^)]*)\\)", "\\1",  tax )
 #   
 #   bips = L[(bipart+2):(length(L)-1)]
 #   blb = L[(blob+2):(bipart-2)]
 #   ntax=length(tax)
 #   tr= rmtree(2,n = ntax,tip.label = 1:ntax)
 #   qTab= quartetTable(tr)
 #   qTab[,(dim(qTab)[2]-3):(dim(qTab)[2]-1)]=qTab[,(dim(qTab)[2]-3):(dim(qTab)[2]-1)]*0 +50
 #   drows=c()
 #   attributes(qTab)$dimnames[[2]] <- c(taxon,'12|34', '13|24', '14|23', '1234') 
 #   
 #   for(i in 1:length(bips))
 #   {
 #     seq=gsub(' ','',bips[i])
 #     bits <- strsplit(seq, "")[[1]]
 #     positions_ones <- which(bits == "1")
 #     positions_zeroes <- which(bits == "0")
 #     
 #     combinationsones <- combn(positions_ones, 2)
 #     combinationszeroes <- combn(positions_zeroes, 2)
 #     
 #     for(ii in 1:(dim(combinationsones)[2]))
 #     {
 #       xx= sort(combinationsones[,ii])
 #       for(jj in 1:(dim(combinationszeroes)[2]))
 #       {
 #         yy=  sort(combinationszeroes[,jj])
 #         target =  c(rep(0,ntax),50,50,50,0)
 #         target[ c(xx,yy)] = 1
 #         drow = which(apply(qTab, 1, function(row) all(row == target)))
 #         drows=c(drows,drow)
 #         if(xx[2]<yy[1] || yy[2]<xx[1]){tag = c(100,25,25,0)}else
 #           if((xx[1]<yy[1] && yy[2]<xx[2]) || (yy[1]<xx[1] && xx[2]<yy[2])) {tag =c(25,25,100,0)}else{tag =  c(25,100,25,0)}
 #         qTab[drow,(dim(qTab)[2]-3):dim(qTab)[2]] = tag
 #       }
 #     }
 #     
 #   }
 #   
 #   
 #   for(i in 1:length( blb))
 #   {
 #     seq=gsub(' ','', blb[i])
 #     bits <- as.numeric(strsplit(seq, "")[[1]])
 #     csize=max(bits)+1
 #     positions_zeroes <- which(bits == 0)
 #     pos=list(positions_zeroes)
 #     for(z in 2:csize)
 #     {
 #       pos[[z]]= which(as.numeric(bits) ==  z-1)
 #     }
 #     
 #     combs= combn(1:csize, 4) 
 #     
 #     for(z in 1:dim(combs)[2])
 #     {
 #       desk = sort(combs[,z])
 #       grd = expand.grid(pos[desk[1]][[1]],pos[desk[2]][[1]], pos[desk[3]][[1]], pos[desk[4]][[1]])
 #       
 #       for(zz in 1:(dim(grd)[1]))
 #       {
 #         
 #         
 #         target =  c(rep(0,ntax),50,50,50,0)
 #         target[ unlist(grd[zz,]) ] = 1
 #         drow = which(apply(qTab, 1, function(row) all(row == target)))
 #         drows=c(drows,drow)
 #         
 #         if(1%in%desk)
 #         {
 #           tag=c(50,50,50,0)
 #           for(zzz in 1:2)
 #           {
 #             xx = sort(unlist(grd[zz,][c(1,2*zzz) ]))
 #             yy = sort(unlist(setdiff(grd[zz,],xx)))
 #             
 #             if(xx[2]<yy[1] || yy[2]<xx[1]){tag = c(25,0,0,0)+tag}else
 #               if((xx[1]<yy[1] && yy[2]<xx[2]) || (yy[1]<xx[1] && xx[2]<yy[2])) {tag =c(0,0,25,0)+tag}else{tag =  c(0,25,0,0)+tag}
 #             
 #           }
 #           qTab[drow,(dim(qTab)[2]-3):dim(qTab)[2]] = tag
 #           
 #           
 #         }else{
 #           
 #           for(zzz in 1)
 #           {
 #             xx = sort(unlist(grd[zz,][c(1,2*zzz) ]))
 #             yy = sort(unlist(setdiff(grd[zz,],xx)))
 #             
 #             if(xx[2]<yy[1] || yy[2]<xx[1]){tag = c(100,25,25,0)}else
 #               if((xx[1]<yy[1] && yy[2]<xx[2]) || (yy[1]<xx[1] && xx[2]<yy[2])) {tag =c(25,25,100,0)}else{tag =  c(25,100,25,0)}
 #             
 #           }
 #           qTab[drow,(dim(qTab)[2]-3):dim(qTab)[2]] = tag
 #         }    
 #         
 #         
 #       }
 #       
 #       
 #     }
 #     
 #     for(sa in 1:length(pos))
 #     {
 #       if(length(pos[[sa]])>1)
 #       {
 #         dd = combn(pos[[sa]],2)
 #         
 #         combined <- unlist(pos)
 #         combined <-setdiff(combined,pos[[sa]])
 #         
 #         dds = combn(combined,2)
 #         
 #         for(isd in 1:(dim(dd)[2]))
 #         {
 #           xx = sort(as.vector(dd[,isd]))
 #           for(isds in 1:(dim(dds)[2]))
 #           {
 #             yy = sort(as.vector(dds[,isds]))
 #             target =  c(rep(0,ntax),50,50,50,0)
 #             target[ c(xx,yy)] = 1
 #             drow = which(apply(qTab, 1, function(row) all(row == target)))
 #             if(length(drow)>0)
 #             {
 #               if(xx[2]<yy[1] || yy[2]<xx[1]){tag = c(100,25,25,0)}else
 #                 if((xx[1]<yy[1] && yy[2]<xx[2]) || (yy[1]<xx[1] && xx[2]<yy[2])) {tag =c(25,25,100,0)}else{tag =  c(25,100,25,0)}
 #               qTab[drow,(dim(qTab)[2]-3):dim(qTab)[2]] = tag
 #             }
 #           }
 #           
 #         }
 #       }
 #       
 #     }
 #     
 #     
 #     
 #   }
 #   
 #   tob =  TINNIK(qTab)
 #   ss=resolveLevel1(tob$ToB,tob$pTable,alpha=0.05,beta=0.95)
 #   return(ss[[1]][[1]])
 # }