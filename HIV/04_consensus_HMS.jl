#!/usr/bin/env julia
#
# Support values for the HMS consensus networks (Huber, Moulton & Spillner
# 2023). Each HMS network is taken as a fixed reference, and support values are
# mapped onto it from the squirrel network set it was built from: its blobs,
# circular orders, hybrid clades and non-redundant bipartitions are each scored
# by how often they occur across those networks. This maps support onto a given
# network -- it is not a bootstrap. Runs for both subtype selections (all, and
# without J/H), one report per run.
#
# The HMS networks are hardcoded below (see HMS_ALL / HMS_NOJH), so this file is
# self-contained -- no separate perl step. Each was produced by 
# HMS-consensus.pl; the exact command is shown above each network. That script
# needs the perl modules Bio::PhyloNetwork and Graph (`cpanm Bio::PhyloNetwork
# Graph`) and is not part of this folder.
#
# Note: HMS networks can contain polytomies (non-binary blobs); PhyloSummaries
# then warns that a blob's circular order is not unique. This is expected and
# does not affect the reported support values.
#
# Requires (versions used): Julia 1.11.7, PhyloSummaries v0.1.1,
# PhyloNetworks v1.3.1, DataFrames v1.8.2. To install:
#   julia> ] add PhyloNetworks DataFrames
#   julia> ] add https://github.com/JuliaPhylo/PhyloSummaries.jl

using PhyloSummaries, PhyloNetworks, DataFrames

# Paths are relative to this script, so the folder can be moved or cloned.
const SAMPLE_DIR = joinpath(@__DIR__, "networks_squirrel")
const OUTPUT_DIR = joinpath(@__DIR__, "networks_consensus")
mkpath(OUTPUT_DIR)

# HMS consensus networks, one per subtype selection. Reproduce each from a
# checkout of the consensus programs, reading the rooted subtype-labelled sample:
#
#   perl HMS-consensus.pl \
#       --in=networks_squirrel/squirrel_rooted_sampling_HIV_all_subtypename.txt \
#       --out=-
const HMS_ALL = "(((K,(F2,F1)I3)I1,(B,D,(A6B)#H1)I6,#H1,A2,H,G,(A3,A1)I4,J,A4,C,A6)I2,N)I5;"

#   perl HMS-consensus.pl \
#       --in=networks_squirrel/squirrel_rooted_sampling_HIV_noJH_subtypename.txt \
#       --out=-
const HMS_NOJH = "(N,((C,(((F2,F1)I7,K)I5,(D,(B,(A6B)#H1)I1)I11)I13)I3,(G,((((A3,A1)I6,(A6,#H1)I12)I4,A4)I2,A2)I14)I9)I10)I8;"

# Each entry is (output name, reference network as newick, squirrel sample). The
# rooted samples are used, as those are the networks the HMS consensus summarises.
const RUNS = [
    ("HMSconsensus_HIV_all",  HMS_ALL,  "squirrel_rooted_sampling_HIV_all_subtypename.txt"),
    ("HMSconsensus_HIV_noJH", HMS_NOJH, "squirrel_rooted_sampling_HIV_noJH_subtypename.txt"),
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

for (name, refnewick, samplefile) in RUNS
    refnet    = readnewick(refnewick)
    netsample = readmultinewick(joinpath(SAMPLE_DIR, samplefile))
    res = blobpartitions_support(netsample, refnet)

    outfile = joinpath(OUTPUT_DIR, "$name.txt")
    open(outfile, "w") do io
        println(io, "Reference: HMS consensus network (hardcoded in this script)")
        println(io, "Sample:    $samplefile")
        println(io, "Networks:  $(length(netsample))")

        println(io, "\n=== Reference network (newick, with internal labels) ===")
        println(io, writenewick(refnet))

        println(io, "\n=== Reference network (newick, without internal labels) ===")
        println(io, writenewick(refnet; internallabel = false))

        # Print every support table the result carries, so nothing is silently
        # dropped if PhyloSummaries adds keys later.
        for key in keys(res)
            print_entry(io, key, res[key])
        end
    end
    println("$name: support from $(length(netsample)) networks -> $outfile")
end
