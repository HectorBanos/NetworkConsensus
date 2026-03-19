
## From Cecile

cd("path/to/directory/with/phylip/files")
files = filter!(x -> occursin(r"\.phy$",x), readdir());
outfile = "full.phy"

ngenes = length(files)
ngenes>0 || error("I didn't find any file in the current folder")
println("found $ngenes files")
phyDict = Dict{String,Vector{String}}() # will map each taxon to its list of sequences
taxa1 = String[] # taxa in the order in which they come in first file
filecounter = 0
rx = r"\s*(\w+)\s+(\S+)\s*$"; # regular expression to extract taxon name & sequence
for file in files
    # global taxa1 # for testing, when run in the REPL
    global filecounter += 1
    open(file) do f # open for reading by default
        readline(f) # read & discard header: "numTaxa numSites"
        counter=0
        for line in eachline(f)
            counter += 1
            m = match(rx, line) # strips away leading & trailing spaces
            m != nothing || continue # to next line. accommodates blank or weird lines
            taxon = m.captures[1]
            if !haskey(phyDict, taxon)
              if filecounter==1
                push!(taxa1, taxon)
                phyDict[taxon] = String[] # starting an empty list of sequences
              else
                @warn "taxon $taxon in alignment $filecounter but not in alignment 1: will be ignored"
              end
            end
            push!(phyDict[taxon], m.captures[2]) # DNA sequence from that file
            # not handled here: taxon in file 1 but not in file i>1
            # counter<2 || break # for testing
        end
    end
end
totallength = unique(sum(length.(v)) for v in values(phyDict));
length(totallength) == 1 ||
  error("concatenated sequences have variable lengths: $totallength")
totallength = totallength[1]

open(outfile, "w") do f
    println(f, length(phyDict), " ", totallength)
    for tax in taxa1 # same order as in file 1
      write(f, tax, " ", string(phyDict[tax]...), "\n")
    end
end
