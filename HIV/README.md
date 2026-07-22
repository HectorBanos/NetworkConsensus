# Consensus methods for phylogenetic networks on HIV-1

Sampling phylogenetic networks from HIV-1 whole-genome sequences with
[squirrel](https://github.com/nholtgrefe/squirrel), summarising them with
consensus methods, and comparing against maximum-likelihood trees from RAxML-NG.

Each analysis is run twice: once on all subtypes, and once with subtypes **J**
and **H** removed.

## Background

The data are the 2023 RIP reference alignment from the LANL HIV Sequence
Database, which provides up to four representative reference sequences per
subtype and CRF. We use the 14 pure group M subtypes — A (subsubtypes A1, A2,
A3, A4, A6), B, C, D, F (F1, F2), G, H, J, K — plus the subtype **N** outgroup:
15 taxa and 53 reference sequences (two to four per taxon), together with the
A6/B CRF (CRF03_AB, labelled `A6B` here).

Because squirrel takes a single aligned sequence per taxon while most taxa have
several references, step 1 runs a resampling analysis of 1,000 replicates to
assess robustness to that choice. Each replicate draws one sequence uniformly at
random from each subtype's references (including the N outgroup), always keeping
the single A6/B recombinant sequence (`03_A6B.RU.97.KAL153_2.AF193276`), then
runs squirrel (`lambda = 0.3`) to infer a level-1
network. Because replicates draw from overlapping subsets the 1,000 networks are
not independent, but together they capture the sampling variability across
reference sequences; steps 2–4 summarise them with consensus methods. Squirrel
is run in two modes — rooted at N, and unrooted (semi-directed, no outgroup).

As a tree-based reference for the non-recombinant signal, step 99 infers
maximum-likelihood trees with RAxML-NG (GTR+GAMMA) on the reference-sequence
alignments (53 taxa, and 46 without J/H), excluding the A6/B CRF.

## Pipeline

Scripts are numbered in the order they run. Step 1 produces the networks the
consensus methods depend on; step 99 is independent of the network side.

| Script | Does | Reads | Writes |
| --- | --- | --- | --- |
| `1_sample_squirrel_networks.py` | Samples 1000 networks per subtype selection, rooted and unrooted | `data/` | `networks_squirrel/` |
| `2_consensus_L1N.jl` | Level-1 consensus network + support | `networks_squirrel/` | `networks_consensus/` |
| `3_consensus_HWC.jl` | Hardwired cluster consensus network (hardcoded) + support | `networks_squirrel/` | `networks_consensus/` |
| `4_consensus_HMS.jl` | HMS consensus network (hardcoded) +  support | `networks_squirrel/` | `networks_consensus/` |
| `99_RaxML-NG.sh` | ML trees + bootstrap support | `data/` | `trees_raxmlng/` |

Steps 2–4 are the three consensus methods being compared. Each writes its
consensus network with support (scored from the squirrel sample) to
`networks_consensus/`, and they are independent of each other. L1N (step 2) is
computed directly in Julia; HWC (step 3) and HMS (step 4) are built externally by
perl and hardcoded in their scripts, each with the command that produced it
shown alongside.

## Folders

- **`data/`** — all alignments. `HIV1_RIP_2023_genome_DNA.fasta` is the source
  squirrel samples from. It also contains sequences
  not used in the analysis. The two `references_*` files are the RAxML alignments —
  subsets of the source that do **not** contain the A6B recombinant.
- **`networks_squirrel/`** — the sampled networks, as eNewick, one per line.
  Each subtype selection is sampled in two rooting modes: rooted at the outgroup
  (`squirrel_rooted_*`) and unrooted/semi-directed with no outgroup
  (`squirrel_unrooted_*`). The `_subtypename` versions (subtype labels) sit in
  this folder; their `_seqname` counterparts (full FASTA taxa names) go in the
  `seqname/` subfolder. The consensus methods read the `_subtypename` files,
  since they require matching leaf labels across replicates.
  `top5_frequent_networks.txt` lists the five most frequent networks per run,
  grouped by isomorphism (rooted labeled isomorphism for the rooted runs,
  semi-directed for the unrooted).
- **`networks_consensus/`** — the consensus networks with support, one
  plain-text report each: the network in Newick (with and without internal
  labels) plus its support tables (blob/partition, hybrid, bipartition, and — for
  HMS/HWC — the full circular-order distribution). `HMSconsensus_*` and
  `HWCconsensus_*` cover the two subtype selections; `L1Nconsensus_*` is the
  level-1 consensus for both selections, each rooted and unrooted — eight files.
- **`trees_raxmlng/`** — `raxml_trees.txt` holds both ML trees with
  bootstrap support and the run details; `logs/<alignment>/` keeps the full
  RAxML-NG output per run.

## Sequence sets

The subtype dictionary `RIP_SUBTYPE_SEQUENCES` in `1_sample_squirrel_networks.py`
is the single source of truth for which sequences are used. It maps each subtype
to its isolates; one isolate per subtype is sampled per replicate.

The two RAxML alignments pair up as:

| Alignment | Taxa | Contents |
| --- | --- | --- |
| `references_all` | 53 | the isolates squirrel samples from, minus the A6B recombinant |
| `references_noJH` | 46 | as above, without J and H |

`A6B` is the circulating recombinant form CRF03_AB, represented by the single
sequence `03_A6B.RU.97.KAL153_2.AF193276`. Unlike the other taxa it is never
resampled — that exact sequence is fixed in every squirrel replicate as the
hybrid taxon. It is excluded from the two `references_*` RAxML alignments, which
are a tree-based analysis.

## Running

Python (step 1):

```bash
pip install physquirrel==1.0.9
python 1_sample_squirrel_networks.py     # slow: 4 x 1000 network inferences
```

Julia (steps 2–4). L1N (`2_consensus_L1N.jl`) computes the level-1 consensus
directly. The HMS and HWC consensus networks are hardcoded in `4_consensus_HMS.jl`
and `3_consensus_HWC.jl` (each with the `perl …-consensus.pl` command that
produced it), so no perl is needed to run them. All three then score their
network against the squirrel sample.

```julia
julia> ] add PhyloNetworks DataFrames
julia> ] add https://github.com/JuliaPhylo/PhyloSummaries.jl
```

```bash
julia 2_consensus_L1N.jl
julia 3_consensus_HWC.jl
julia 4_consensus_HMS.jl
```

Regenerating the HMS/HWC networks themselves needs 
`HMS-consensus.pl` / `HWC-consensus.pl` (perl modules `Bio::PhyloNetwork` and
`Graph`); the exact commands are in the two script headers.

RAxML-NG (step 99), expected at `~/miniforge3/bin/raxml-ng`:

```bash
./99_RaxML-NG.sh                         # slow: 2 ML searches with bootstrapping
```

## Notes

- **Outgroup.** Subtype `N` roots the rooted squirrel networks and the consensus;
  the unrooted networks are sampled with no outgroup (squirrel returns a
  semi-directed network, serialised to eNewick with an arbitrary rooting). RAxML-NG
  infers unrooted trees (GTR+G is time-reversible, so its `--outgroup` is a
  drawing option only); root them on `N` afterwards if needed.
- **Versions used.** Python 3.11 · physquirrel 1.0.9 · Julia 1.11.7, PhyloSummaries 0.1.1,
  PhyloNetworks 1.3.1, DataFrames 1.8.2 · RAxML-NG 2.0.2 (conda/miniforge3) ·
  phylox 1.1.2, phylozoo 0.2.6 (network counting). 
