library(MSCquartets)
library(SiPhyNetwork)
path_to_Julia = '/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia'
base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"
location <- file.path(base_dir, "sim_out")
auxfile_location <- file.path(base_dir, "auxfiles")
setwd(auxfile_location)

lambda=1; mu=0.2; nu=0.2;
inheritance_fxn <- make.beta.draw(2, 2)
gen_dist_fxn = make.stepwise(probs = c(1,0),distances = c(1.1,Inf))
hybridproportions = c(0.25,0.55,0.2)
networks=c()
scale=0.75
k=1
while(k<21)
{
      # net<-sim.bdh.taxa.ssa( n=20, numbsim=1, frac = .6,
                       # lambda=lambda, mu=mu, nu=nu, hybprops=hybridproportions,
                       # hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = inheritance_fxn, complete =T)[[1]]
         net<-sim.bdh.taxa.gsa(n=20,m=21,numbsim=1,
                         lambda=lambda, mu=mu, nu=nu, hybprops=hybridproportions,
                          hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = inheritance_fxn, complete = T)[[1]]
      plot(net)
         # net<-sim.bdh.age(age = 10, numbsim=1,
      #                   lambda=lambda, mu=mu, nu=nu, hybprops=hybridproportions,
      #                   hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = inheritance_fxn ) 
 
      if(!is.numeric(net) )
          {
        if( length(net$reticulation)>0)
         {
         
          lvl <- getNetworkLevel(net)
       
           net$edge.length = scale*net$edge.length/mean(net$edge.length)
           netnew=write.net(net)
           setwd(auxfile_location)
           
           J = readLines('mastertob.jl')
           J[2] = paste0('cd(','\"',location,'\"',')')
           J[3] = paste0('net = ', '\"',netnew,'\"')
           write(J, file = 'tob.jl')
           cmd = paste0(path_to_Julia,' ',location,'/tob.jl')
           system(cmd)
          # edgelabels(round(net$edge.length,2), bg="black", col="white", font=2)
           tob=read.tree(file="topo.tree")
           # level=as.numeric(readLines("lvl1.txt")[1])
          # plot(tob)
          deg=ape::degree(tob)
          if(lvl==1 && ((sum(deg[deg[,1]>3,2])>3) ||  (sum(deg[,1]>5)>=2)) ) 
          {
            mean(net$edge.length)
            netnew=write.net(net)
            networks[k]=netnew
            plot(k,k)
            k=k+1
            cat( paste0(netnew,"\n"), file = "lvl1networks_newrun.txt", append = TRUE)
            cat( paste0(readLines('topo.tree'),"\n"), file = "lvl1networks_tob_newrun.txt", append = TRUE)
          }
        }
      }
}
   