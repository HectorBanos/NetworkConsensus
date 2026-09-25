#!/usr/bin/env julia
#
# Where the root sits in the two rooted squirrel samples. Writes
# networks_squirrel/root_frequencies.txt (pass another directory as the first
# argument to write elsewhere).
#
# In each rooted network the outgroup N is a child of the root; the root's other
# child -- the attached node -- induces a tripartition of the taxa via its three
# incident edges: the outgroup side {N}, and the two taxon sets reachable below
# its two child edges (hardwired clusters, PhyloNetworks' hardwiredcluster).
# Networks are grouped by this tripartition and the TOP most frequent are listed.
#
# Requires (version used): Julia 1.11.7, PhyloNetworks v1.3.1.

using PhyloNetworks, Printf

const DIR = joinpath(@__DIR__, "networks_squirrel")
const OUTFILE = joinpath(length(ARGS) >= 1 ? ARGS[1] : DIR, "root_frequencies.txt")
const OUTGROUP = "N"
const TOP = 5
const RUNS = [("Rooted, all subtypes",    "squirrel_rooted_sampling_HIV_all_subtypename.txt"),
              ("Rooted, without J and H", "squirrel_rooted_sampling_HIV_noJH_subtypename.txt")]

# "N | set | set" for one network, the two sets in alphabetical order.
function tripartition(net)
    taxa = sort(tiplabels(net))
    kids = getchildren(getroot(net))
    out = findfirst(k -> k.leaf && k.name == OUTGROUP, kids)
    att = kids[3 - out]                       # the attached node
    sets = sort([taxa[hardwiredcluster(e, taxa)] for e in att.edge if getparent(e) === att])
    return join([OUTGROUP; join.(sets, ",")], " | ")
end

open(OUTFILE, "w") do io
    println(io, """
Where the root sits in the two rooted squirrel samples.

In each rooted network the outgroup N is a child of the root; the root's other
child -- the attached node -- induces a tripartition of the taxa via its three
incident edges: the outgroup side {N}, and the two taxon sets reachable below
its two child edges. We list the five most frequent tripartitions per run,
written  N | set | set.  If the attached node is a blob node (the reticulation
cycle runs through it), the hybrid taxon is reachable down both child edges and
so appears in both sets. Reachable sets are hardwired clusters, computed with
PhyloNetworks. The hybrid is whichever taxon squirrel inferred as reticulate in
that replicate -- usually A6B, occasionally another (e.g. C).""")
    for (title, file) in RUNS
        lines = filter(!isempty ∘ strip, readlines(joinpath(DIR, file)))
        n = length(lines)
        counts = Dict{String,Int}(); order = String[]     # order = first-seen, for ties
        for l in lines
            t = tripartition(readnewick(l))
            haskey(counts, t) || push!(order, t)
            counts[t] = get(counts, t, 0) + 1
        end
        top = first(sort(order; by = t -> -counts[t], alg = MergeSort), TOP)
        s = sum(counts[t] for t in top)
        println(io, "\n", "="^72, "\n", title, "\n", "="^72)
        println(io, "Input:    ", file, "\nNetworks: ", n, "\nDistinct: ", length(order))
        @printf(io, "Top %d:    %d/%d (%.1f%%)\n", TOP, s, n, 100s / n)
        for (i, t) in enumerate(top)
            @printf(io, "\n#%d  %d/%d  (%.1f%%)\n%s\n", i, counts[t], n, 100counts[t] / n, t)
        end
    end
end
println("wrote ", OUTFILE)
