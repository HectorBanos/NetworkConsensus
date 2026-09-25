#!/usr/bin/env julia
#
# All clades that occur in the rooted squirrel network samples (all subtypes,
# and without J/H). A clade is counted in a network if it is a hardwired
# cluster of that network (PhyloNetworks' hardwiredclusters); we list every
# non-trivial clade found in at least one network, most frequent first. Writes
# networks_squirrel/clade_frequencies.txt (pass another directory as the first
# argument to write elsewhere).
#
# Trivial clades are left out: a single taxon, everything but one taxon (e.g.
# all taxa except the outgroup), and the full set.
#
# Requires (version used): Julia 1.11.7, PhyloNetworks v1.3.1.

using PhyloNetworks, Printf

const DIR = joinpath(@__DIR__, "networks_squirrel")
const OUTFILE = joinpath(length(ARGS) >= 1 ? ARGS[1] : DIR, "clade_frequencies.txt")
const RUNS = [("Rooted, all subtypes",    "squirrel_rooted_sampling_HIV_all_subtypename.txt"),
              ("Rooted, without J and H", "squirrel_rooted_sampling_HIV_noJH_subtypename.txt")]

# Non-trivial clades of one network, each as a sorted vector of taxa.
function clades(net)
    taxa = sort(tiplabels(net))
    M = hardwiredclusters(net, taxa)          # one row per internal edge, 0/1 per taxon
    out = Set{Vector{String}}()
    for i in 1:size(M, 1)
        c = taxa[M[i, 2:end-1] .== 1]
        2 <= length(c) <= length(taxa) - 2 && push!(out, c)
    end
    return out
end

open(OUTFILE, "w") do io
    println(io, """
Frequencies of all clades occurring in the rooted squirrel network samples.

For each sampling run, the count is how many of the input networks contain the
split "clade | rest" -- i.e. have that clade as a hardwired cluster of the
network (computed with PhyloNetworks). Trivial clades (a single taxon, all but
one taxon, the full taxon set) are left out.""")
    for (title, file) in RUNS
        lines = filter(!isempty ∘ strip, readlines(joinpath(DIR, file)))
        n = length(lines)
        counts = Dict{Vector{String},Int}()
        for l in lines, c in clades(readnewick(l))
            counts[c] = get(counts, c, 0) + 1
        end
        freq = sort([(join(c, ","), k) for (c, k) in counts]; by = x -> (-x[2], x[1]))
        w = maximum(length(x[1]) for x in freq)
        println(io, "\n", "="^72, "\n", title, "\n", "="^72)
        println(io, "Input:    ", file, "\nNetworks: ", n, "\nClades:   ", length(freq), "\n")
        for (c, k) in freq
            @printf(io, "%s | rest : %5d/%d  (%5.1f%%)\n", rpad(c, w), k, n, 100k / n)
        end
    end
end
println("wrote ", OUTFILE)
