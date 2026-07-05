using PhyloNetworks
using PhyloCoalSimulations
 
cd("/Users/006512736/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents")
num_gts = 10000;
newick = "(C:1.0,(((B:1.0,(A:1.0)#H1:1.0::0.5):1.0,((F:1.0,#H1:1.0::0.5):1.0,(E:1.0,(D:1.0)#H2:1.0::0.5):1.0):1.0):1.0,#H2:1.0::0.5):1.0);"
graph_name = "test"
ntwk = readTopology(newick);
gt_fname1 = string("sampled_gtrees_",graph_name)
gts1 = simulatecoalescent(ntwk, num_gts, 1; inheritancecorrelation=0.0 );
writeMultiTopology(gts1,gt_fname1)
exit()
