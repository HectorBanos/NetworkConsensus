#!/usr/bin/env julia
#
# Level-1 consensus networks from the squirrel network samples.
# Runs on both rooting modes -- rooted (rooted networks, outgroup N) and
# unrooted (unrooted/semi-directed networks, no outgroup) -- for both subtype
# selections (all subtypes, and without J/H), writing one report per run.
#
# Requires (versions used): Julia 1.11.7, PhyloSummaries v0.1.1,
# PhyloNetworks v1.3.1, DataFrames v1.8.2. To install:
#   julia> ] add PhyloNetworks DataFrames
#   julia> ] add https://github.com/JuliaPhylo/PhyloSummaries.jl

using PhyloSummaries, PhyloNetworks, DataFrames

# Paths are relative to this script, so the folder can be moved or cloned.
const INPUT_DIR  = joinpath(@__DIR__, "networks_squirrel")
const OUTPUT_DIR = joinpath(@__DIR__, "networks_consensus")
mkpath(OUTPUT_DIR)

# Consensus needs matching leaf labels across replicates, hence the
# subtype-labelled networks. The rooted runs use "N" as the outgroup (as in the
# squirrel sampling); the unrooted runs pass no outgroup (outgroup = nothing),
# matching the unrooted network sets. Each entry is (output name, input, outgroup).
const RUNS = [
    ("L1Nconsensus_rooted_HIV_all",    "squirrel_rooted_sampling_HIV_all_subtypename.txt",    "N"),
    ("L1Nconsensus_rooted_HIV_noJH",   "squirrel_rooted_sampling_HIV_noJH_subtypename.txt",   "N"),
    ("L1Nconsensus_unrooted_HIV_all",  "squirrel_unrooted_sampling_HIV_all_subtypename.txt",  nothing),
    ("L1Nconsensus_unrooted_HIV_noJH", "squirrel_unrooted_sampling_HIV_noJH_subtypename.txt", nothing),
]

# Print a result entry: tables in full (no row/column truncation), plain
# vectors (e.g. :taxa) as a list.
function print_entry(io, title, val)
    println(io, "\n=== $title ===")
    if val isa AbstractVector{<:AbstractString}
        println(io, join(val, ", "))
        return
    end
    df = DataFrame(val, copycols = false)
    show(io, df; allrows = true, allcols = true, truncate = 0)
    println(io)
    println(io, "Shape: $(nrow(df)) rows × $(ncol(df)) columns")
end

for (name, infile, outgroup) in RUNS
    netsample = readmultinewick(joinpath(INPUT_DIR, infile))
    res = consensus_level1network(netsample; outgroup = outgroup, suppressinfo = true)

    outfile = joinpath(OUTPUT_DIR, "$name.txt")
    open(outfile, "w") do io
        println(io, "Input:    $infile")
        println(io, "Networks: $(length(netsample))")
        println(io, "Outgroup: $(isnothing(outgroup) ? "none (unrooted)" : outgroup)")

        println(io, "\n=== Consensus network (newick, with internal labels) ===")
        println(io, writenewick(res[:net]))

        println(io, "\n=== Consensus network (newick, without internal labels) ===")
        println(io, writenewick(res[:net]; internallabel = false))

        # Print every table the result carries, so nothing is silently dropped
        # if PhyloSummaries adds keys later.
        for key in keys(res)
            key === :net && continue
            print_entry(io, key, res[key])
        end
    end
    println("$name: consensus of $(length(netsample)) networks -> $outfile")
end
