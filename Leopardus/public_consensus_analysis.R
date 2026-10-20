# R code for reproducing bootstrap analyses in 
#       `Consensus methods for trees of blobs and level-1 phylogenetic networks` 
#        by K. Wicke, H. Banos, N. Holtgrefe, E. S. Allman, C. Ane, and O. Gascuel
#
# September 2026

# This code is a template for creating bootstrap replicates, both trees of blobs
# and networks for downstream tree of blobs and L1N construction.  Either
# TINNiK or ECToBlob can be used tree of blob construction, but TINNiK is used in demo.
#
# The user must set a working directory below.  In this directory a file called 
# `gt_file.newick' # stores the original gene trees in extended Newick format.  
# Also below, the user # should assign an integer to a variable `runNo`.   
# In this example runNo = 1.  When # runing this script, directories will be 
# created where the bootstrap output is stored.
#
# Specifially,
#
#   run1/bs_tob/         :  TINNiK_tob_bs has bootstrap ToBs
#   run1/bs_net/         :  multiplicities_bs_nets_TINNiK     vector of multiplicities
#                           net_bs_TINNiK                     bootstap networks using TINNiK
#                           ntwks_with_weights_TINNiK.csv     csv file for PhyloSummaries input (networks and weights)   
#   run1/bs_vec/         :  Leopardus_r2                      data_name_rX where X is number of replicates in bootstrap analysis
#   run1/ref_ToB_ntwk/   :  ToB, level-1 network files.  Plus file with settings.

rm(list=ls())

library(MSCquartets)
library(tictoc)

# **** change to your working directory *****
main_dir = "~/Desktop/github_files"
setwd(main_dir)

# **** your gene tree file is assumed to be in the working
#      directory with the name gt_file.newick                ****
gt_fname = "gt_file.newick"

runNo = 1

# to suppress messages for bootstrap analysis, set to true
quiet = TRUE

# set the number of bootstrap replicates
num_bs_reps = 31

# set the runNo directory
run_dir = file.path(main_dir,paste0("run",runNo))

# create directories if needed
{
  if (!dir.exists(run_dir))  dir.create(run_dir)
  tob_dir = file.path(run_dir,"bs_tob")
  tob_vec_dir = file.path(run_dir,"bs_vec")
  net_dir = file.path(run_dir,"bs_net")
  ref_dir = file.path(run_dir,"ref_ToB_ntwk")
  
  if (!dir.exists(tob_dir))  dir.create(tob_dir)
  if (!dir.exists(tob_vec_dir))  dir.create(tob_vec_dir)
  if (!dir.exists(net_dir))  dir.create(net_dir)
  if (!dir.exists(ref_dir))  dir.create(ref_dir)
}

# Boolean control variables
doWriteRefTob_and_ntwk = TRUE   # compute the TINNiK ToB and NANUQp level-1 network and store to files
doBootstrapToB = TRUE           # compute the bootstrap ToBs and store to files
do_NANUQp = TRUE                # compute the NANUQp resolutions of bootstrap ToBs

writeToFile = TRUE

# set the values for the hypothesis tests used by TINNiK, ECToBlob, and NANUQ+
alpha_value=5e-29
beta_value=0.95

# In the paper the Leopardus gene trees  of Lescroart et al, [29] were used.  These
# are publicly available.  If taxon names are not supplied, then they are selected
# from the first gene tree.  This can be a problem if there are taxa missing on this
# first tree.  For safety, you can set the taxa names here.

# this data_name is used in filenames for storing bootstrap data
# data_name = "Your_choice"

# Example in paper
data_name = "Leopardus"
# Comment or uncomment as appropriate.
tNames = c("Lcolocola_chile","Lcolocola_3" ,"Lguigna_chile",
           "Ltigrinus_pardinoides_32454","Ltigrinus_pardinoides_32451","Ltigrinus_oncilla_13",
           "Lguttulus_brazil","Ltigrinus_emiliae_699","Ltigrinus_emiliae_6",
           "Lgeoffroyi_3","Lpardalis_6","Lpardalis_brazil",
           "Lpardalis_texas","Ljacobita_chile" ,"Lwiedii_brazil",
           "Puma_concolor" )


tic(" Total elapsed time.")

tic("Read gene trees")
gts=read.tree(file= gt_fname)
toc()

num_gts = length(gts)
cat(paste0("Number of gts is: ",num_gts,"\n"))

# Get the reference ToB and NANUQp network(s) for the gene trees in gt_file.newick
if (doWriteRefTob_and_ntwk) {
  
  if (writeToFile) {
    if (!dir.exists(ref_dir))  dir.create(ref_dir)
  }
  
  # do TINNiK.  An alternative is to use ECToBlob to get the tree of blobs.
  {
    # get TINNiK tree
    tic("Computing TINNiK tree of blobs for gene trees")
    output = TINNIK(gts, alpha = alpha_value, beta = beta_value, plot=FALSE, taxonnames = tNames)
    pTable = output$pTable
    TINNiK_ToB = output$ToB
    toc()
  }
  
  # now get reference level-1 network(s)
  tic("Computing NANUQp level-1 networks using tree of blobs for gene trees")
  N_TINNiK = resolveLevel1(TINNiK_ToB, pTable, alpha = alpha_value, beta = beta_value, plot = 0, fullResMax = 20)
  toc()
  
  # networks from NANUQp using TINNiK ToB
  nets_TINNiK <- N_TINNiK[[1]]
  num_TINNiK_nets <- length(nets_TINNiK)
  
  cat(paste0("\nNumber of level-1 networks returned by NANUQp with TINNiK tree of blobs as input: ",num_TINNiK_nets,"\nTest levels are alpha: ", alpha_value,
             " and beta: ", beta_value,"\n"))
  
  # for writing the ToB
  if (writeToFile) {
    ToB_fname <- file.path(ref_dir,paste0("TINNiK_ref_ToB.tre"))
    file.create(ToB_fname)
    write.tree(TINNiK_ToB,file=ToB_fname)
    
    # write out a settings file for easy reference 
    settings_fname <- file.path(ref_dir,paste0("settings_",data_name,".txt"))
    cat(paste0(data_name,":  ", as.character(Sys.time())), "\n", file = settings_fname)
    cat("alpha: ", alpha_value,"\n", file = settings_fname, append = TRUE)
    cat("beta: ", beta_value,"\n", file = settings_fname, append = TRUE)
    cat("num taxa: ",Ntip(TINNiK_ToB),"\n", file = settings_fname, append = TRUE)
    
    # now write out networks
    nets_fname <- file.path(ref_dir,paste0("TINNiK_ref_ntwk.newick"))
    fh <- file(nets_fname, "w")
    if (!is.null(nets_TINNiK)) {
      for (row in nets_TINNiK) {
        writeLines(paste(row, collapse = " "), fh)
      }
    } else {
      writeLines(write.tree(TINNiK_ToB),fh)
      print("No NANUQp nets with TINNiK")
    }
    close(fh)
    
  }
}

# Bootstrap gene trees to get bootstrap ToBs.  Then save to files.  
# If do_NANUQp is true, then also do bootstrap for level-1 NANUQp networks, 
if (doBootstrapToB) {
  
  # ECToBlob is a second method for getting ToB, arguably better than TINNiK
  # method_names = <- c("TINNiK", "EB")
  method_names <- c("TINNiK")
  
  # setwd, then return to main directory once finished
  setwd(run_dir)
  
  TINNiK_ToB_fname = file.path(tob_dir,"TINNiK_tob_bs")
  
  # set up storage
  if (do_NANUQp) {
    
    # Since NANUQp can return ties, prepare data structures to keep track of them
    num_opt_ntwks <- vector("list", length(method_names))
    names(num_opt_ntwks) <- method_names
    
    # for the networks
    opt_ntwks <- vector("list", length(method_names))
    names(opt_ntwks) <- method_names
    
    # fnames for storing the networks
    net_fname_TINNiK <- file.path(net_dir, "net_bs_TINNiK")
    
    # fname for storing the multiplicities
    mult_vec_fname = file.path(net_dir,"multiplicities_bs_nets_TINNiK")
  }
  
  cat(paste0("\nStarting to compute bootstrap samples, ToBs, level-1 resolutions for ",num_bs_reps," replicates.\n"))
  tic("Total time for bootstrap replicates")
  # cat("\n")
  
  # get file names to be of same length
  for (rep in 1:num_bs_reps) {
    
    # for storing the resampling vectors
    fname = file.path(tob_vec_dir, sprintf("%s_r%s", data_name, num_bs_reps))
    
    # sample gene trees with replacement 
    cols = sample.int(num_gts, num_gts, replace = TRUE)
    
    write(cols, file = fname, ncolumns = length(cols), append = (rep > 1))
    
    # bootstrap gene tree sample
    boot_gts = gts[cols]
    
    # do TINNiK first
    {    
      
      if (quiet) {
        suppressMessages({output=TINNIK(boot_gts,alpha=alpha_value,beta=beta_value,plot=FALSE,taxonnames=tNames)})
      } else {
        output=TINNIK(boot_gts,alpha=alpha_value,beta=beta_value,plot=FALSE,taxonnames=tNames)
      }
      
      bootstrap_pTable=output$pTable
      ToB = output$ToB
      
      # label internal nodes
      TINNiK_ToB = labelIntNodes(ToB, plot=FALSE)
      
      # write TINNiK_ToBs with internal nodes labeled to file
      if (rep==1) {
        cat(write.tree(TINNiK_ToB), "\n", file = TINNiK_ToB_fname)
      }
      if (rep>1) {
        cat(write.tree(TINNiK_ToB), "\n", file = TINNiK_ToB_fname, append = TRUE)
      } 
      
    }
    
    # Get the NANUQp bootstrap level-1 networks.  
    if (do_NANUQp) {
      
      if (rep < 10) {
        str = paste0("00",rep)
      }
      if ((rep > 9) && (rep < 100)) {
        str = paste0("0",rep)
      }
      if (rep > 99) {
        str = rep
      }
      
      # TINNiK example
      {
        m_name <- "TINNiK"
        
        ToB = TINNiK_ToB
        
        # Only do resolveLevel1 if it is not binary, and
        # if ToB is binary, write it to the file directly; 
        if (is.binary(ToB)) {
          
          n_newick = write.tree(ToB, file = "")
          numTiedRes = 1 # no tied resolutions
          
          num_opt_ntwks[[m_name]][rep] <- numTiedRes
          opt_ntwks[[m_name]] <- c(opt_ntwks[[m_name]],n_newick)
          
          if (rep==1) {
            write(n_newick,file=net_fname_TINNiK)
          }
          if (rep>1) {
            write(n_newick,file=net_fname_TINNiK,append = TRUE)
          } 
        } else {
          if (quiet) {
            suppressMessages({N = resolveLevel1(ToB, bootstrap_pTable, alpha = alpha_value, beta = beta_value, plot = 0, fullResMax = 20)})
          } else {
            N = resolveLevel1(ToB, bootstrap_pTable, alpha = alpha_value, beta = beta_value, plot = 0, fullResMax = 20)
          }
          
          # check for tied resolutions and store rep number, number of ties, 
          # and eventually write all optimal ntwks to file
          numTiedRes = length(N[[1]])
          num_opt_ntwks[[m_name]][rep] <- numTiedRes
          opt_ntwks[[m_name]] <- c(opt_ntwks[[m_name]],unlist(N[[1]]))
          
        }
        
        write(
          unlist(N[[1]]),
          file = net_fname_TINNiK,
          append = (rep > 1)
        )
        
      } 
    } # end TINNiK resolveLevel1
    
    # could do same with ECToBlob
    
    if ((rep %% 10)==0) cat("Finished replicate",rep,"\n")
  } # end rep loop
  
  
  
  # after all the replicates are finished, write out multiplicities of network to file
  if (do_NANUQp) {
    
    cat("\nWriting out vector of multiplicities to file.\n")
    
    # only TINNiK in this template
    m_name = "TINNiK"
    
    # write out vector with number of optimal networks to file
    write(num_opt_ntwks[["TINNiK"]], ncolumns = length(cols), file = mult_vec_fname)
    
    # now write out file for use in Julia's PhyloSummaries
    fname <- file.path(net_dir,paste0("ntwks_with_weights_",m_name,".csv"))
    cat("net,weight\n", file = fname)
    
    nets <- opt_ntwks[[m_name]]
    mult <- num_opt_ntwks[[m_name]]
    weights <- rep(ifelse(mult == 1, "1", paste0("1/", mult)), mult)
    for (i in seq_along(nets)) {
      cat(sprintf('"%s",%s\n', nets[i], weights[i]),file = fname, append = TRUE)
    }
    rm(list=c("mult","nets"))
  }
  
  toc()
  cat("\n")
} 



# return to main directory once finished
setwd(main_dir)

toc()
