#######################################################
########## Create subsequences of length 300 ##########
#######################################################


#For all netwokrs each gene we obtain a subsequence creating a shorter MSA

library(ape)   

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"
location <- file.path(base_dir, "sim_out")
setwd(location)

lens <- read.csv("seq_lens.csv", row.names = 1)$x    

starts <- cumsum(c(1, lens[-length(lens)]))           
ends   <- cumsum(lens)   
n_sites <-c()
seq_lens<-c()
for(xx in 1:20)
{
  gts_filetag = paste0('RealisticSim_400gts_Net',xx)
  
  aln <- read.dna(paste0("concatalign", gts_filetag, ".fasta"), format = "fasta", as.matrix = TRUE)
  
  
  set.seed(42)    
  for(i in 1:400)
  {
    n_sites[i] <- ceiling(lens[i]/3)
    seq_lens[i] =  ceiling((lens/3)[i])
  }
  
  sampled_cols <- lapply(seq_along(lens), function(i) {
    gene_cols <- starts[i]:ends[i]
    
    # sample(gene_cols, n_sites[i], replace = FALSE)
     head(gene_cols, n_sites[i]) ##this is we want the first n_sites
  })
  
  
  all_cols <- sort(unlist(sampled_cols))
  
  
  aln_sub <- aln[, all_cols]
  
  
  write.dna(aln_sub, file = paste0("concatalign", gts_filetag, "_third.fasta"),
            format = "fasta", nbcol = -1, colsep = "")
}



write.csv(seq_lens, file=paste0(location,'/seq_lens_third.csv'))


