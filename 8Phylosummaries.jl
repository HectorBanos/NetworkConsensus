using PhyloSummaries
using PhyloNetworks
nets=String["nets"] 
base_dir = expanduser("~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2")
sim_dir = joinpath(base_dir, "sim_out")
cd(sim_dir)

nets = []
for I in 1:20
    netfile = "rooted_network_newick_concatalignRealisticSim_400gts_Net$(I).txt"
    netsample = readmultinewick(netfile)
    netsamp = netsample[2:101]
    res_l1 = consensus_level1network(netsamp)
    con = res_l1[:net]
    push!(nets, writenewick(con, internallabel=false))
end


filename = "Consensus_outgroup_full.txt"
open(filename, "w") do io
    for line in nets
        println(io, line)
    end
end


nets = []
for I in 1:20
    netfile = "rooted_network_newick_concatalignRealisticSim_400gts_Net$(I)_third.txt"
    netsample = readmultinewick(netfile)
    netsamp = netsample[2:101]
    res_l1 = consensus_level1network(netsamp)
    con = res_l1[:net]
    push!(nets, writenewick(con, internallabel=false))
end


filename = "Consensus_outgroup_third.txt"
open(filename, "w") do io
    for line in nets
        println(io, line)
    end
end




nets = []
for I in 1:20
    netfile = "unrooted_network_newick_concatalignRealisticSim_400gts_Net$(I).txt"
    netsample = readmultinewick(netfile)
    netsamp = netsample[2:101]
    res_l1 = consensus_level1network(netsamp)
    con = res_l1[:net]
    push!(nets, writenewick(con, internallabel=false))
end


filename = "Consensus_nooutgroup_full.txt"
open(filename, "w") do io
    for line in nets
        println(io, line)
    end
end


nets = []
for I in 1:20
    netfile = "unrooted_network_newick_concatalignRealisticSim_400gts_Net$(I)_third.txt"
    netsample = readmultinewick(netfile)
    netsamp = netsample[2:101]
    res_l1 = consensus_level1network(netsamp)
    con = res_l1[:net]
    push!(nets, writenewick(con, internallabel=false))
end


filename = "Consensus_nooutgroup_third.txt"
open(filename, "w") do io
    for line in nets
        println(io, line)
    end
end

 