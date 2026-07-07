using PhyloNetworks
cd("/Users/006512736/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/simulation")
net = 
net1 = readnewick(net)
ttop = treeofblobs(net1)
writeTopology(ttop,"topo.tree")
exit()
