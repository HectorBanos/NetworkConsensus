 # Repository Overview

This repository contains all the materials required to replicate the analyses presented in *Consensus methods for trees of blobs and level-1 phylogenetic networks*, including the simulation pipeline, simulation analyses, and empirical data analyses (HIV and *Leopardus*).

The repository consists of five main directories and 16 files (numbered from **0** to **15**).

## Directory Structure

### 1. `HIV/`

Contains all files required to reproduce the analyses of the HIV dataset.

### 2. `Leopardus/`

Contains all files required to reproduce the analyses of the *Leopardus* dataset.

### 3. `sim_out/`

Contains all simulated datasets, inferred networks, and consensus networks used in the simulation study.

For each of the 20 simulated networks, the directory contains:

* Gene tree samples of 400 genes, for example:

  * `sampled_gtrees_RealisticSim_400gts_NetX`
  * where `X` denotes the network number.

* Concatenated multiple sequence alignments (MSAs):

  * `concatalignRealisticSim_400gts_NetX.fasta`

* Sequence lengths for the 400 simulated genes:

  * `seq_lens.csv`

* Reduced MSAs, obtained by retaining only one-third of each gene:

  * `concatalignRealisticSim_400gts_NetX_third.fasta`

* Sequence lengths for the reduced alignments:

  * `seq_lens_third.csv`
    
### Replicate Subfolders

All network inference and consensus results are organised into **10 independent replicates**. Each replicate corresponds to an independent run of Squirrel on fresh random seeds, allowing variability in the inference to be assessed. The replicates are numbered 1 through 10.

For each replicate `i` (where `i` ranges from 1 to 10), four subfolders are created inside `sim_out/`:

| Subfolder | Rooting | Alignment |
|---|---|---|
| `rootedfull`*`i`* | Rooted (outgroup provided) | Full alignment |
| `rootedthird`*`i`* | Rooted (outgroup provided) | One-third alignment |
| `unrootedfull`*`i`* | Unrooted (no outgroup) | Full alignment |
| `unrootedthird`*`i`* | Unrooted (no outgroup) | One-third alignment |

For example, `rootedfull3/` contains all outputs from the third replicate of the rooted, full-alignment analysis, and `rootedthird7/` contains all outputs from the seventh replicate of the rooted, one-third-alignment analysis.

Each subfolder contains the following files for all 20 simulated networks:

* The Squirrel output files (101 networks each; see below).
* The consensus networks produced by PhyloSummaries.jl (one per setting; see below).
* The HMS and HWC consensus networks (one file per simulated network; see below).
* The processed results used for figures and tables (`.Rdata` files produced by scripts 9–12).
* The summary spreadsheets (`.xlsx` files produced by script 13).

### Squirrel Network Inference

Each replicate subfolder contains the Squirrel output for all 20 simulated networks in that setting. For example, in `rootedfull`*`i`*`/`:

* `rooted_network_newick_concatalignRealisticSim_400gts_NetX.txt`

contains 101 networks:

* The first network is the Squirrel network inferred from the original dataset;
* The remaining 100 networks correspond to bootstrap replicates.

In the unrooted subfolders the filenames begin with `unrooted_` instead of `rooted_`, and in the one-third subfolders the filenames include the suffix `_third` before `.txt`.

In total, each replicate subfolder contains 20 Squirrel output files (one per simulated network), and there are 40 subfolders across 10 replicates and 4 analysis settings, for a grand total of 800 Squirrel output files.

### Consensus Networks

Each replicate subfolder contains one consensus file produced by PhyloSummaries.jl for the corresponding analysis setting:

* `rootedfull`*`i`*`/` → `Consensus_outgroup_full.txt`
* `rootedthird`*`i`*`/` → `Consensus_outgroup_third.txt`
* `unrootedfull`*`i`*`/` → `Consensus_nooutgroup_full.txt`
* `unrootedthird`*`i`*`/` → `Consensus_nooutgroup_third.txt`

Each file contains 20 Newick networks, one per line. Line *i* corresponds to the consensus network for simulated network *i*.

### HMS and HWC Consensus Networks

Each replicate subfolder contains the HMS and HWC consensus networks for all 20 simulated networks in that setting. For example, in `rootedfull`*`i`*`/`:

* `HMS_outgroup_full_NetX.txt`
* `HWC_outgroup_full_NetX.txt`

and in `rootedthird`*`i`*`/`:

* `HMS_outgroup_third_NetX.txt`
* `HWC_outgroup_third_NetX.txt`

where `X` denotes the simulated network number (1–20).

---

### 4. `auxfiles/`

Contains auxiliary files required to generate the simulated datasets.

Although many files are used by the simulation pipeline, the most relevant are:

* `HMS-consensus.pl`: Perl implementation of the HMS consensus method.
* `HWC-consensus.pl`: Perl implementation of the HWC consensus method.
* `XXX.lisp`: Lisp implementation of the consensus method introduced in this work.

---

### 5. `DataforEstimates/`

Contains all files used to estimate the simulation parameters.

These include:

* the sequence alignments for the 10 primate genes,
* the Paioinini dataset,
* the IQ-TREE output files,
* an R script used to process the IQ-TREE output and estimate the parameters used for the simulations.

## Repository Files

The repository also contains 16 numbered files (`0`–`15`), which implement the simulation pipeline and reproduce all analyses presented in the paper.

For all scripts, the only modifications required before execution are to update the file paths so that they match your local installation. Depending on the script, this may include:

* the path to the cloned repository,
* the path to the IQ-TREE executable,
* the path to the Julia executable.

### Simulation Pipeline (Files 0–8)

#### `0Sim_Networks.R`

Generates the 20 simulated networks satisfying the conditions described in the paper (i.e., networks containing either at least two cycles of size greater than four or at least four cycles of size greater than three).

This is the most computationally intensive step of the pipeline. Since all simulated networks are already provided in the repository, we do **not** recommend rerunning this script unless new simulations are desired.

#### `1Sim_seq.R`

Simulates, for each network:

* 400 gene trees,
* the corresponding sequence alignments for each gene, and
* the concatenated multiple sequence alignment (MSA).

#### `2Subset_seq.R`

Produces the reduced MSAs by retaining only one-third of the sites from each simulated gene.

#### `3SQUIRREL_rooted_full.py`

Runs Squirrel on the full concatenated alignments while providing the root, producing both the inferred network and 100 bootstrap networks. This script is run once per replicate, writing its outputs to the corresponding `rootedfull`*`i`*`/` subfolder.

#### `4SQUIRREL_rooted_partial.py`

Runs Squirrel on the one-third alignments while providing the root. Outputs are written to `rootedthird`*`i`*`/`.

#### `5SQUIRREL_unrooted_full.py`

Runs Squirrel on the full concatenated alignments without providing the root. Outputs are written to `unrootedfull`*`i`*`/`.

#### `6SQUIRREL_unrooted_partial.py`

Runs Squirrel on the one-third alignments without providing the root. Outputs are written to `unrootedthird`*`i`*`/`.

Scripts 3–6 are each executed 10 times (once per replicate). An R wrapper script is provided that automates the 10 repetitions for all four scripts, creating the replicate subfolders, copying the required MSA and sequence-length files into each one, and invoking the appropriate Python script.

#### `7Huber_HWC_cmdlines.txt`

Contains the terminal commands required to compute the HMS and HWC consensus networks for all simulation settings. These commands are run for each replicate subfolder.

#### `8PhyloSummaries.jl`

Julia implementation used to compute the consensus networks introduced in this work from the Squirrel bootstrap networks. This script is run for each of the 40 replicate subfolders (10 replicates × 4 analysis settings), reading the Squirrel bootstrap networks from the subfolder and writing the consensus network file to the same subfolder.

### Analysis and Figure Generation (Files 9–15)

Scripts 9–13 are each run across all 10 replicates, reading inputs from and writing outputs to the replicate-specific subfolders. Scripts 14–15 aggregate the results across all 10 replicates, averaging the metrics before generating the final figures.

* `9Process_edges_outgroup.R` – Processes edge-based metrics for the rooted analyses, looping over all 10 replicates and writing results to `rootedfull`*`i`*`/` and `rootedthird`*`i`*`/`.
* `10Process_edges_nooutgroup.R` – Processes edge-based metrics for the unrooted analyses, writing results to `unrootedfull`*`i`*`/` and `unrootedthird`*`i`*`/`.
* `11Process_blobs_outgroup.R` – Processes blob, hybrid, and circular order-based metrics for the rooted analyses, across all 10 replicates.
* `12Process_blobs_nooutgroup.R` – Processes blob, hybrid, and circular order-based metrics for the unrooted analyses, across all 10 replicates.
* `13CreateExcels.R` – Compiles the processed results into summary spreadsheets, writing one pair of `.xlsx` files into each of the 40 replicate subfolders.
* `14CreatePlots.R` – Aggregates results across all 10 replicates, computes averages, and generates the main figures presented in the paper.
* `15CreatePlots_disttoinput.R` – Aggregates results across all 10 replicates, computes averages, and generates figures summarizing the distance of the inferred consensus networks to the true simulated networks.
