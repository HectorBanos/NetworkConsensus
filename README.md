# Consensus Methods for Trees of Blobs and Level-1 Phylogenetic Networks

This repository contains all materials needed to replicate the analyses in *Consensus methods for trees of blobs and level-1 phylogenetic networks*, including:

- the simulation pipeline,
- the simulation analyses, and
- the empirical analyses (HIV and *Leopardus*).

 
---

## Repository Structure

```
.
├── HIV/                     # HIV empirical analysis
├── Leopardus/               # Leopardus empirical analysis
├── sim_out/                 # Simulated data, inferred and consensus networks
├── networksjl_adjusted      # The 20 model networks used in the simulation study
├── HMS-consensus.pl         # HMS consensus method (Perl)
├── HWC-consensus.pl         # HWC consensus method (Perl)
└── L1Nconsensus_Lisp.lsp                 # L1N Consensus method introduced in this work (Lisp)
 
```

---
## Consensus Method Implementations

| File | Description |
|---|---|
| `HMS-consensus.pl` | Perl implementation of the HMS consensus method |
| `HWC-consensus.pl` | Perl implementation of the HWC consensus method |
| `L1Nconsensus_Lisp.lsp` | Lisp implementation of the consensus method introduced in this work |

---
## Empirical Data

| Directory | Contents |
|---|---|
| `HIV/` | All files required to reproduce the analyses of the HIV dataset. |
| `Leopardus/` | All files required to reproduce the analyses of the *Leopardus* dataset. |

---

## Model Networks

**`networksjl_adjusted`** contains the **20 model (true) networks** used in the simulation study. Gene trees and alignments in `sim_out/` were simulated from these networks, and the `X`-th network in this file corresponds to `NetX` throughout the rest of the repository.

---
 
## Simulation Output (`sim_out/`)

Contains all simulated datasets, inferred networks and consensus networks used in the simulation study. Throughout, `X` denotes the simulated network number (1–20).

### Simulated Data

For each of the 20 model networks:

| File | Description |
|---|---|
| `sampled_gtrees_RealisticSim_400gts_NetX` | Sample of 400 simulated gene trees |
| `concatalignRealisticSim_400gts_NetX.fasta` | Concatenated multiple sequence alignment (MSA) |
| `seq_lens.csv` | Sequence lengths of the 400 simulated genes |
| `concatalignRealisticSim_400gts_NetX_third.fasta` | Reduced MSA, keeping one-third of each gene |
| `seq_lens_third.csv` | Sequence lengths for the reduced alignments |

### Replicate Subfolders

All network inference and consensus results are organised into **10 independent replicates** (`i = 1, …, 10`). Each replicate is an independent bootstrap sample to run of Squirrel, allowing the variability of inference to be assessed.

For each replicate `i`, there are four subfolders:

| Subfolder | Rooting | Alignment |
|---|---|---|
| `rootedfull<i>/` | Rooted (outgroup provided) | Full |
| `rootedthird<i>/` | Rooted (outgroup provided) | One-third |
| `unrootedfull<i>/` | Unrooted (no outgroup) | Full |
| `unrootedthird<i>/` | Unrooted (no outgroup) | One-third |

*Example:* `rootedfull3/` holds the third replicate of the rooted, full-alignment analysis; `rootedthird7/` holds the seventh replicate of the rooted, one-third-alignment analysis.

Each subfolder contains, for all 20 simulated networks:

- Squirrel output files (101 networks each)
- PhyloSummaries.jl consensus networks (one file per setting)
- HMS and HWC consensus networks (one file per simulated network)

### Squirrel Network Inference

File naming pattern:

```
<rooting>_network_newick_concatalignRealisticSim_400gts_NetX[_third].txt
```

- `<rooting>` is `rooted` or `unrooted`
- `_third` is present only in the one-third subfolders

Example: `rooted_network_newick_concatalignRealisticSim_400gts_NetX.txt`

Each file contains **101 networks**:

1. Line 1 — the Squirrel network inferred from the original data
2. Lines 2–101 — networks inferred from 100 bootstrap replicates

Each subfolder has 20 such files (one per simulated network): **40 subfolders × 20 files = 800 Squirrel output files** in total.

### PhyloSummaries.jl Consensus Networks

| Subfolder | Consensus file |
|---|---|
| `rootedfull<i>/` | `Consensus_outgroup_full.txt` |
| `rootedthird<i>/` | `Consensus_outgroup_third.txt` |
| `unrootedfull<i>/` | `Consensus_nooutgroup_full.txt` |
| `unrootedthird<i>/` | `Consensus_nooutgroup_third.txt` |

Each file contains 20 Newick networks, one per line; **line `X` is the consensus network for simulated network `X`**.

### HMS and HWC Consensus Networks

One file per simulated network and method:

| Subfolder | Files |
|---|---|
| `rootedfull<i>/` | `HMS_outgroup_full_NetX.txt`, `HWC_outgroup_full_NetX.txt` |
| `rootedthird<i>/` | `HMS_outgroup_third_NetX.txt`, `HWC_outgroup_third_NetX.txt` |

---


 
