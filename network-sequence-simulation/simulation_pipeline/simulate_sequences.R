##################################################
########## Sequence simulation pipeline ##########
##################################################
library(rlist)
library(phytools)
library(MSCquartets)
library(MCMCpack)
library(ape)
library(beepr)
library(igraph)
library(phangorn)
library(pvclust)
library(BAMMtools)
library(castor)
library(ape)
library(phangorn)
library(adephylo)
library(TreeTools)

# Set working directory and directories for outut and auxiliary files. 
base_dir <- "path/to/your/project"
location <- file.path(base_dir, "simulation_output")
auxfile_location <- file.path(base_dir, "simulation_auxfile")

# Set paths to Julia and IQ-TREE on your system
path_to_Julia  <- "/path/to/julia"
path_to_IQtree <- "/path/to/iqtree3"

# Number of gene trees
num_gts = 400

# Random seed
set.seed(1961)

# Create random gene lengths centered around 600 bp and save lengths to a file.
seq_lens = round(300+300*rgamma(num_gts, 5,5),0)
write.csv(seq_lens, file=paste0(location,'/seq_lens.csv'))

# Read in a file of phylogenetic networks.
networks = readLines("networksjl_adjusted.txt")

# Read in realistic tree diameter values (based on Papionini).
dtreespapio=as.numeric(read.csv(file="papionini_diameters.csv")[[1]])

# Loop over all input networks.
for(xx in 1:length(networks))
{
        # Choose a filetag to name output files.
        gts_filetag = paste0('RealisticSim_400gts_outgroup_Net',xx)
        
        # Create a Julia script that simulates random gene trees for the current
        # network using PhyloCoalSimulations.jl.
        setwd(auxfile_location)
        JAne = readLines('ConcatANEte.jl')
        
        gts_filetagge = paste0('Gtrees_Net',xx)
        J = readLines('masterphylocoal.jl')
        J[4] = paste0('cd(','\"',location,'\"',')')
        J[5] = paste0('num_gts = ',num_gts,';')
        J[6] = paste0('newick = ','\"',networks[xx],'\"')
        J[7] = paste0('graph_name =','\"',gts_filetag,'\"')
      
        # Run the resulting Julia script.
        setwd(location)
        write(J, file = 'PhyCoalSim.jl')
        cmd = paste0(path_to_Julia,' ',location,'/PhyCoalSim.jl')
        system(cmd)

        # Read in the random gene trees generated in the previous step.
        gts=read.tree(file=paste0("sampled_gtrees_",gts_filetag))
        
        # Set seed.
        set.seed(1961)
        
        # Sample a diameter from the collection of realistic (Papionini) 
        # diameters.
        randdiam = sample( dtreespapio ,num_gts,replace = T) 
        
        # Set seed.
        set.seed(1961)
        
        seque = paste0(auxfile_location,'/papiofull_out.phylip')
        
        # Set seed based on the current network number. 
        set.seed(xx)
        
        # Create a set of random integer seeds to be used when calling IQ-TREE
        # (Alisim) below. Each gene tree will have its own seed.
        seeds=ceiling(2000*runif(num_gts))
 
        # To break ultrametricity in the gene trees, each edge length is scaled
        # by an independent draw from a log-normal distribution with mean 1 and
        # standard deviation 0.1278835. This standard deviation was selected so
        # that, on average, the ratio of the maximum to minimum root-to-tip
        # distances across taxa is approximately 1.2.
        sigma = 0.1278835 
        
        maxratiosim=c()
        diamsim = c()
        
        # Loop over gene trees.
        for(i in 1:num_gts)
        {
          # Break ultrametricity as described above.
          maxratiosim[i]=5
          while( maxratiosim[i]>4.9)
          {
                ultramebreak = rlnorm(length(gts[[i]]$edge.length),-0.5 * sigma^2, sigma)
                gts[[i]]$edge.length = gts[[i]]$edge.length*ultramebreak
                rescale = randdiam [i] * 1/max(cophenetic(gts[[i]]))
                gts[[i]]$edge.length=gts[[i]]$edge.length*rescale
               
                diamsim[i]=max(cophenetic( gts[[i]]))
                distances_all_castor <- get_all_distances_to_root(gts[[i]])[1:length(gts[[i]]$tip.label)]
                maxratiosim[i] = max(distances_all_castor)/min(distances_all_castor)
          }
          
          # Add an outgroup to the tree.
          new_tree <- AddTip(tree = gts[[i]], where = 0, label = "U", lengthBelow = diamsim[i]/2)
          new_tree$edge.length[new_tree$edge.length==0] = 2*diamsim[i]/2
          gts[[i]]<- new_tree
          
          # Prepare IQ-TREE/Alisin call to simulate sequence data under a GTR
          # model with rates {1.02895, 4.07930, 0.62855, 1.32375, 3.74745}, base
          # frequencies {0.28250, 0.21125, 0.20285, 0.29685}, gamma shape = 2.95.
          write.tree(gts[[i]],file=paste0("ssampled_gtrees_",i,'.treefile'))
          cmdalisim=paste0(path_to_IQtree,' ','--alisim sseq',i,' -seed ',seeds[i],' -redo -m GTR{1.02895/4.07930/0.62855/1.32375/3.74745}+F{0.28250/0.21125/0.20285/0.29685}+G{2.95} -te ',paste0(location,"/ssampled_gtrees_",i,'.treefile'),' --seqtype  DNA  --length ',seq_lens[i],' -blfix -af phy   -s ',seque )
          
          # Call IQ-TREE.
          iqali=system(cmdalisim,intern = T);
          
          # Remove auxiliary files.
          file.remove(paste0("ssampled_gtrees_",i,'.treefile'))
          
          # Show progress.
          if(i==100 || i==200 || i==300){print(paste0("Net ",xx," Tree ",i))}
    }
        
        # Create a Julia script to concatenate the sequence data for the 
        # individual gene trees. 
        JAne[4] = paste0('cd(','\"',location,'\"',')')
        setwd(location)
        write(JAne, file = 'ConcatANEte_sloc.jl')
        
        # Run the Julia script.
        cmd = paste0(path_to_Julia,' ',location,'/ConcatANEte_sloc.jl')
        system(cmd)
        
        files <- paste0("sseq", 1:num_gts, ".phy")
        file.remove(files)
        
        # Create file for concatenated data.
        file.rename('full.phy', paste0('concatalign',gts_filetag,'.txt'))
        paste0('concatalign',gts_filetag,'.txt')
        
        # Remove auxiliary files.
        file.remove( file = 'PhyCoalSim.jl')
        file.remove( file = 'ConcatANEte_sloc.jl')
       
}