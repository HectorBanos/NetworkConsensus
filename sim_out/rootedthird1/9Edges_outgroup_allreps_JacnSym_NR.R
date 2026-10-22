#########################################################################
#################### Process edge info   rooted case ####################
#########################################################################

### Computes BOTH Jaccard and symmetric (RF) distances in one pass,
### saving separate .Rdata files for each.

library(ape)
library(SiPhyNetwork)
library(grid)
library(gridExtra)
library(TreeDist)
library(MSCquartets)
library('Quartet')

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"

source(paste0(base_dir,"/auxfiles/Finalcode_bootsupport.R"))

## ---------------------------------------------------------------
## Helper: splits restricted to edges with degree-3 on BOTH sides
## ---------------------------------------------------------------
degree3_split_names <- function(tree) {
  
  M <- split_matrix(tree)[-1, , drop = FALSE]
  if (nrow(M) == 0) return(character(0))
  
  deg     <- table(c(tree$edge[, 1], tree$edge[, 2]))
  n_tips  <- length(tree$tip.label)
  n_nodes <- max(tree$edge)
  
  descendants <- vector("list", n_nodes)
  for (i in seq_len(n_tips)) descendants[[i]] <- i
  for (i in rev(seq_len(nrow(tree$edge)))) {
    parent <- tree$edge[i, 1]; child <- tree$edge[i, 2]
    descendants[[parent]] <- unique(c(descendants[[parent]], descendants[[child]]))
  }
  
  keep <- sapply(seq_len(nrow(M)), function(r) {
    row   <- M[r, ]; target <- row[1]
    side1 <- which(row == target)
    matched_edge <- NA_integer_
    for (e in seq_len(nrow(tree$edge))) {
      desc <- descendants[[tree$edge[e, 2]]]
      if (setequal(desc, side1) || setequal(desc, setdiff(seq_len(n_tips), side1))) {
        matched_edge <- e; break
      }
    }
    if (is.na(matched_edge)) return(FALSE)
    deg_parent <- as.integer(deg[as.character(tree$edge[matched_edge, 1])])
    deg_child  <- as.integer(deg[as.character(tree$edge[matched_edge, 2])])
    isTRUE(deg_parent == 3) && isTRUE(deg_child == 3)
  })
  
  names_all <- sapply(seq_len(nrow(M)), function(i) {
    row <- M[i, ]; target <- row[1]
    paste(names(row[row == target]), collapse = ",")
  })
  names_all[keep]
}

## ---------------------------------------------------------------
## Comparison functions — all use degree3_split_names
## ---------------------------------------------------------------

## Returns c(a, b, c): matched, only-in-tree1, only-in-tree2
split_counts <- function(names1, names2) {
  a <- length(intersect(names1, names2))
  b <- length(setdiff(names1, names2))
  c <- length(setdiff(names2, names1))
  c(a, b, c)
}

edge_dist <- function(tree, boot_tree) {
  split_counts(degree3_split_names(tree), degree3_split_names(boot_tree))
}

jaccard_dist <- function(names1, names2) {
  a <- length(intersect(names1, names2))
  b <- length(setdiff(names1, names2))
  c <- length(setdiff(names2, names1))
  union_size <- a + b + c
  if (union_size == 0) return(0)
  (b + c) / union_size
}

symmetric_dist <- function(names1, names2) {
  a <- length(intersect(names1, names2))
  b <- length(setdiff(names1, names2))
  c <- length(setdiff(names2, names1))
  denom <- b + c          # = |names1| + |names2|  when a subtracted
  if (denom == 0) return(0)
  (b + c) / (length(names1) + length(names2))   # normalised RF
}

edge_supp <- function(tree, boot_trees) {
  names1 <- degree3_split_names(tree)
  edges  <- rep(0, length(names1))
  for (j in seq_along(names1)) {
    nam <- names1[j]
    for (i in seq_along(boot_trees)) {
      names2   <- degree3_split_names(boot_trees[[i]])
      edges[j] <- (nam %in% names2) + edges[j]
    }
  }
  names1[which(edges > 50)]
}

## RF.DIFM variants — return c(distance, a, b, c)
RF.DIFM.jaccard <- function(boot_tree, tree, boot_trees) {
  names2 <- edge_supp(tree, boot_trees)
  names1 <- degree3_split_names(boot_tree)
  a <- length(intersect(names1, names2))
  b <- length(setdiff(names1, names2))
  c <- length(setdiff(names2, names1))
  union_size <- a + b + c
  jac <- if (union_size == 0) 0 else (b + c) / union_size
  c(jac, a, b, c)
}

RF.DIFM.symmetric <- function(boot_tree, tree, boot_trees) {
  names2 <- edge_supp(tree, boot_trees)
  names1 <- degree3_split_names(boot_tree)
  a <- length(intersect(names1, names2))
  b <- length(setdiff(names1, names2))
  c <- length(setdiff(names2, names1))
  denom <- length(names1) + length(names2)
  sym <- if (denom == 0) 0 else (b + c) / denom
  c(sym, a, b, c)
}

## ---------------------------------------------------------------
## Main loop — one pass, both distances computed simultaneously
## ---------------------------------------------------------------
setwd(paste0(base_dir,"/auxfiles"))

networks = readLines("networksjl_outgroup_adjusted_labelled.txt")
loopsize = 20
tobs     = readLines('lvl1networks_tob_letters.txt')

n_reps <- 1

for (rep_i in 1:n_reps) {
  
  message("=== Edges outgroup: replicate ", rep_i, " ===")
  
  ## ---- rootedfull<rep_i> ----
  location <- file.path(base_dir, "sim_out", paste0("rootedfull", rep_i))
  setwd(location)
  
  MM_jac <- matrix(rep(0, 32*20), nrow = 20)
  MM_sym <- matrix(rep(0, 32*20), nrow = 20)
  
  for (xx in 1:loopsize) {
    
    TrueNet <- read.evonet(text = networks[xx])
    TrueToB <- treeOfBlobs(TrueNet)
    TrueToB$edge.length <- rep(1, length(TrueToB$edge.length))
    if (xx < 10) { xp <- paste0('', xx) } else { xp <- xx }
    
    sLvl1 <- read.evonets(file = paste0(location, '/Consensus_outgroup_full.txt'))[[xx]]
    
    flag  <- 1
    hLvl1 <- readLines(paste0(location, '/Huber_outgroup_full_Net', xx, '.txt'))
    if (str_detect(hLvl1, "#")) { hLvl1 <- read.evonet(text = hLvl1) } else { hLvl1 <- read.tree(text = hLvl1); flag <- 0 }
    
    flag1  <- 1
    hwLvl1 <- readLines(paste0(location, '/HWC_outgroup_full_Net', xx, '.txt'))
    if (str_detect(hwLvl1, "#")) { hwLvl1 <- read.evonet(text = hwLvl1) } else { hwLvl1 <- read.tree(text = hwLvl1); flag1 <- 0 }
    
    sqLvl1 <- read.nets(file = paste0(location, '/rooted_network_newick_concatalignRealisticSim_400gts_Net', xp, '.txt'))[1][[1]]
    
    sToB  <- treeOfBlobs(sLvl1)
    if (flag  == 1) { hToB  <- treeOfBlobs(hLvl1)  } else { hToB  <- unroot(hLvl1)  }
    if (flag1 == 1) { hwToB <- treeOfBlobs(hwLvl1) } else { hwToB <- unroot(hwLvl1) }
    sqToB <- treeOfBlobs(sqLvl1)
    sToB$edge.length <- rep(1, length(sToB$edge.length))
    
    sboot_sampled_lvl1 <- read.nets(file = paste0(location, '/rooted_network_newick_concatalignRealisticSim_400gts_Net', xp, '.txt'))[2:101]
    sboot_sampled_tob  <- list()
    for (i in 1:100) {
      flag_b <- 1
      if (is.null(sboot_sampled_lvl1[[i]]$reticulation)) { flag_b <- 0 }
      if (flag_b == 1) { sboot_sampled_tob[[i]] <- treeOfBlobs(sboot_sampled_lvl1[[i]]) } else { sboot_sampled_tob[[i]] <- sboot_sampled_lvl1[[i]] }
    }
    
    ## Pre-compute split names used repeatedly
    sn_true  <- degree3_split_names(TrueToB)
    sn_s     <- degree3_split_names(sToB)
    sn_sq    <- degree3_split_names(sqToB)
    sn_h     <- degree3_split_names(hToB)
    sn_hw    <- degree3_split_names(hwToB)
    sn_boot  <- lapply(sboot_sampled_tob, degree3_split_names)
    sn_supp  <- edge_supp(sToB, sboot_sampled_tob)   # majority-supported splits
    
    ## --- branch counts (unfiltered) ---
    numbran   <- dim(TrueToB$edge)[1] - 21
    numbrans  <- dim(sToB$edge)[1]    - 21
    numbransq <- dim(sqToB$edge)[1]   - 21
    numbranh  <- dim(hToB$edge)[1]    - 21
    numbranhw <- dim(hwToB$edge)[1]   - 21
    numbranm  <- length(sn_supp)
    
    ## --- edge-set distances ---
    sed  <- split_counts(sn_true, sn_s)
    sqed <- split_counts(sn_true, sn_sq)
    hed  <- split_counts(sn_true, sn_h)
    hwed <- split_counts(sn_true, sn_hw)
    
    ## majority RF counts (a,b,c) — same for both distance types
    sn_supp_names <- sn_supp
    a_m <- length(intersect(sn_true, sn_supp_names))
    b_m <- length(setdiff(sn_true, sn_supp_names))
    c_m <- length(setdiff(sn_supp_names, sn_true))
    med <- c(a_m, b_m, c_m)
    
    ## --- Jaccard distances ---
    srf_j   <- jaccard_dist(sn_true, sn_s)
    sqrf_j  <- jaccard_dist(sn_true, sn_sq)
    hrf_j   <- jaccard_dist(sn_true, sn_h)
    hwrf_j  <- jaccard_dist(sn_true, sn_hw)
    union_m <- a_m + b_m + c_m
    mrf_j   <- if (union_m == 0) 0 else (b_m + c_m) / union_m
    
    ## --- Symmetric (normalised RF) distances ---
    srf_s   <- symmetric_dist(sn_true, sn_s)
    sqrf_s  <- symmetric_dist(sn_true, sn_sq)
    hrf_s   <- symmetric_dist(sn_true, sn_h)
    hwrf_s  <- symmetric_dist(sn_true, sn_hw)
    denom_m <- length(sn_true) + length(sn_supp_names)
    mrf_s   <- if (denom_m == 0) 0 else (b_m + c_m) / denom_m
    
    ## --- Bootstrap loop (compute both distances simultaneously) ---
    srfbs_j  <- 0; sqrfbs_j <- 0; hrfbs_j <- 0; hwrfbs_j <- 0; mrfbs_j <- 0; Trfbs_j <- 0
    srfbs_s  <- 0; sqrfbs_s <- 0; hrfbs_s <- 0; hwrfbs_s <- 0; mrfbs_s <- 0; Trfbs_s <- 0
    
    for (i in 1:100) {
      sn_bi <- sn_boot[[i]]
      
      ## majority RF for bootstrap sample i
      a_bi <- length(intersect(sn_bi, sn_supp_names))
      b_bi <- length(setdiff(sn_bi, sn_supp_names))
      c_bi <- length(setdiff(sn_supp_names, sn_bi))
      u_bi <- a_bi + b_bi + c_bi
      d_bi <- length(sn_bi) + length(sn_supp_names)
      
      srfbs_j  <- jaccard_dist(sn_bi, sn_s)   / 100 + srfbs_j
      sqrfbs_j <- jaccard_dist(sn_bi, sn_sq)  / 100 + sqrfbs_j
      hrfbs_j  <- jaccard_dist(sn_bi, sn_h)   / 100 + hrfbs_j
      hwrfbs_j <- jaccard_dist(sn_bi, sn_hw)  / 100 + hwrfbs_j
      mrfbs_j  <- (if (u_bi == 0) 0 else (b_bi + c_bi) / u_bi) / 100 + mrfbs_j
      Trfbs_j  <- jaccard_dist(sn_bi, sn_true)/ 100 + Trfbs_j
      
      srfbs_s  <- symmetric_dist(sn_bi, sn_s)   / 100 + srfbs_s
      sqrfbs_s <- symmetric_dist(sn_bi, sn_sq)  / 100 + sqrfbs_s
      hrfbs_s  <- symmetric_dist(sn_bi, sn_h)   / 100 + hrfbs_s
      hwrfbs_s <- symmetric_dist(sn_bi, sn_hw)  / 100 + hwrfbs_s
      mrfbs_s  <- (if (d_bi == 0) 0 else (b_bi + c_bi) / d_bi) / 100 + mrfbs_s
      Trfbs_s  <- symmetric_dist(sn_bi, sn_true)/ 100 + Trfbs_s
    }
    
    row_common <- c(numbran, numbrans, numbransq, numbranh, numbranhw, numbranm,
                    sed, sqed, hed, hwed, med)
    
    MM_jac[xx, ] <- c(row_common, srf_j, sqrf_j, hrf_j, hwrf_j, mrf_j,
                      Trfbs_j, srfbs_j, sqrfbs_j, hrfbs_j, hwrfbs_j, mrfbs_j)
    MM_sym[xx, ] <- c(row_common, srf_s, sqrf_s, hrf_s, hwrf_s, mrf_s,
                      Trfbs_s, srfbs_s, sqrfbs_s, hrfbs_s, hwrfbs_s, mrfbs_s)
  }
  
  MM1 <- MM_jac; save(MM1, file = paste0(location, '/EdgeAnalysis_full_outgroup_jaccard.Rdata'))
  MM1 <- MM_sym; save(MM1, file = paste0(location, '/EdgeAnalysis_full_outgroup_symmetric.Rdata'))
  
  ## ---- rootedthird<rep_i> ----
  location <- file.path(base_dir, "sim_out", paste0("rootedthird", rep_i))
  setwd(location)
  
  MM_jac <- matrix(rep(0, 32*20), nrow = 20)
  MM_sym <- matrix(rep(0, 32*20), nrow = 20)
  
  for (xx in 1:loopsize) {
    
    TrueNet <- read.evonet(text = networks[xx])
    TrueToB <- treeOfBlobs(TrueNet)
    TrueToB$edge.length <- rep(1, length(TrueToB$edge.length))
    if (xx < 10) { xp <- paste0('', xx) } else { xp <- xx }
    
    flags <- 1
    sLvl1 <- readLines(paste0(location, '/Consensus_outgroup_third.txt'))[xx]
    if (str_detect(sLvl1, "#")) { sLvl1 <- read.evonet(text = sLvl1) } else { sLvl1 <- read.tree(text = sLvl1); flags <- 0 }
    
    flag  <- 1
    hLvl1 <- readLines(paste0(location, '/Huber_outgroup_third_Net', xx, '.txt'))
    if (str_detect(hLvl1, "#")) { hLvl1 <- read.evonet(text = hLvl1) } else { hLvl1 <- read.tree(text = hLvl1); flag <- 0 }
    
    flag1  <- 1
    hwLvl1 <- readLines(paste0(location, '/HWC_outgroup_third_Net', xx, '.txt'))
    if (str_detect(hwLvl1, "#")) { hwLvl1 <- read.evonet(text = hwLvl1) } else { hwLvl1 <- read.tree(text = hwLvl1); flag1 <- 0 }
    
    sqLvl1 <- read.nets(file = paste0(location, '/rooted_network_newick_concatalignRealisticSim_400gts_Net', xp, '_third.txt'))[1][[1]]
    
    if (flags == 1) { sToB  <- treeOfBlobs(sLvl1)  } else { sToB  <- unroot(sLvl1)  }
    if (flag  == 1) { hToB  <- treeOfBlobs(hLvl1)  } else { hToB  <- unroot(hLvl1)  }
    if (flag1 == 1) { hwToB <- treeOfBlobs(hwLvl1) } else { hwToB <- unroot(hwLvl1) }
    sqToB <- treeOfBlobs(sqLvl1)
    sToB$edge.length <- rep(1, length(sToB$edge.length))
    
    sboot_sampled_lvl1 <- read.nets(file = paste0(location, '/rooted_network_newick_concatalignRealisticSim_400gts_Net', xp, '_third.txt'))[2:101]
    sboot_sampled_tob  <- list()
    for (i in 1:100) {
      flag_b <- 1
      if (is.null(sboot_sampled_lvl1[[i]]$reticulation)) { flag_b <- 0 }
      if (flag_b == 1) { sboot_sampled_tob[[i]] <- treeOfBlobs(sboot_sampled_lvl1[[i]]) } else { sboot_sampled_tob[[i]] <- sboot_sampled_lvl1[[i]] }
    }
    
    sn_true  <- degree3_split_names(TrueToB)
    sn_s     <- degree3_split_names(sToB)
    sn_sq    <- degree3_split_names(sqToB)
    sn_h     <- degree3_split_names(hToB)
    sn_hw    <- degree3_split_names(hwToB)
    sn_boot  <- lapply(sboot_sampled_tob, degree3_split_names)
    sn_supp  <- edge_supp(sToB, sboot_sampled_tob)
    
    numbran   <- dim(TrueToB$edge)[1] - 21
    numbrans  <- dim(sToB$edge)[1]    - 21
    numbransq <- dim(sqToB$edge)[1]   - 21
    numbranh  <- dim(hToB$edge)[1]    - 21
    numbranhw <- dim(hwToB$edge)[1]   - 21
    numbranm  <- length(sn_supp)
    
    sed  <- split_counts(sn_true, sn_s)
    sqed <- split_counts(sn_true, sn_sq)
    hed  <- split_counts(sn_true, sn_h)
    hwed <- split_counts(sn_true, sn_hw)
    
    sn_supp_names <- sn_supp
    a_m <- length(intersect(sn_true, sn_supp_names))
    b_m <- length(setdiff(sn_true,   sn_supp_names))
    c_m <- length(setdiff(sn_supp_names, sn_true))
    med <- c(a_m, b_m, c_m)
    
    srf_j   <- jaccard_dist(sn_true, sn_s)
    sqrf_j  <- jaccard_dist(sn_true, sn_sq)
    hrf_j   <- jaccard_dist(sn_true, sn_h)
    hwrf_j  <- jaccard_dist(sn_true, sn_hw)
    union_m <- a_m + b_m + c_m
    mrf_j   <- if (union_m == 0) 0 else (b_m + c_m) / union_m
    
    srf_s   <- symmetric_dist(sn_true, sn_s)
    sqrf_s  <- symmetric_dist(sn_true, sn_sq)
    hrf_s   <- symmetric_dist(sn_true, sn_h)
    hwrf_s  <- symmetric_dist(sn_true, sn_hw)
    denom_m <- length(sn_true) + length(sn_supp_names)
    mrf_s   <- if (denom_m == 0) 0 else (b_m + c_m) / denom_m
    
    srfbs_j  <- 0; sqrfbs_j <- 0; hrfbs_j <- 0; hwrfbs_j <- 0; mrfbs_j <- 0; Trfbs_j <- 0
    srfbs_s  <- 0; sqrfbs_s <- 0; hrfbs_s <- 0; hwrfbs_s <- 0; mrfbs_s <- 0; Trfbs_s <- 0
    
    for (i in 1:100) {
      sn_bi <- sn_boot[[i]]
      a_bi <- length(intersect(sn_bi, sn_supp_names))
      b_bi <- length(setdiff(sn_bi, sn_supp_names))
      c_bi <- length(setdiff(sn_supp_names, sn_bi))
      u_bi <- a_bi + b_bi + c_bi
      d_bi <- length(sn_bi) + length(sn_supp_names)
      
      srfbs_j  <- jaccard_dist(sn_bi, sn_s)    / 100 + srfbs_j
      sqrfbs_j <- jaccard_dist(sn_bi, sn_sq)   / 100 + sqrfbs_j
      hrfbs_j  <- jaccard_dist(sn_bi, sn_h)    / 100 + hrfbs_j
      hwrfbs_j <- jaccard_dist(sn_bi, sn_hw)   / 100 + hwrfbs_j
      mrfbs_j  <- (if (u_bi == 0) 0 else (b_bi + c_bi) / u_bi) / 100 + mrfbs_j
      Trfbs_j  <- jaccard_dist(sn_bi, sn_true) / 100 + Trfbs_j
      
      srfbs_s  <- symmetric_dist(sn_bi, sn_s)    / 100 + srfbs_s
      sqrfbs_s <- symmetric_dist(sn_bi, sn_sq)   / 100 + sqrfbs_s
      hrfbs_s  <- symmetric_dist(sn_bi, sn_h)    / 100 + hrfbs_s
      hwrfbs_s <- symmetric_dist(sn_bi, sn_hw)   / 100 + hwrfbs_s
      mrfbs_s  <- (if (d_bi == 0) 0 else (b_bi + c_bi) / d_bi) / 100 + mrfbs_s
      Trfbs_s  <- symmetric_dist(sn_bi, sn_true) / 100 + Trfbs_s
    }
    
    row_common <- c(numbran, numbrans, numbransq, numbranh, numbranhw, numbranm,
                    sed, sqed, hed, hwed, med)
    
    MM_jac[xx, ] <- c(row_common, srf_j, sqrf_j, hrf_j, hwrf_j, mrf_j,
                      Trfbs_j, srfbs_j, sqrfbs_j, hrfbs_j, hwrfbs_j, mrfbs_j)
    MM_sym[xx, ] <- c(row_common, srf_s, sqrf_s, hrf_s, hwrf_s, mrf_s,
                      Trfbs_s, srfbs_s, sqrfbs_s, hrfbs_s, hwrfbs_s, mrfbs_s)
  }
  
  MM2 <- MM_jac; save(MM2, file = paste0(location, '/EdgeAnalysis_third_outgroup_jaccard.Rdata'))
  MM2 <- MM_sym; save(MM2, file = paste0(location, '/EdgeAnalysis_third_outgroup_symmetric.Rdata'))
  
  message("Completed replicate ", rep_i, " of ", n_reps)
}

message("All ", n_reps, " replicates complete (edges, outgroup).")