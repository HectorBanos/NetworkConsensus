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

### Squirrel Network Inference

For each MSA, the directory contains the network inferred by Squirrel together with the **100 bootstrap networks.

For example,

* `rooted_network_newick_concatalignRealisticSim_400gts_NetX.txt`

contains 101 networks:

* The first network is the Squirrel network inferred from the original dataset;
* The remaining 100 networks correspond to bootstrap replicates.

Two inference settings are included:

* Rooted: the root was provided to Squirrel (filenames beginning with `rooted`).
* Unrooted: no root was provided to Squirrel (filenames beginning with `unrooted`).

Both settings are available for:

* the full alignments, and
* The third alignments.

In total, there are 80 Squirrel output files:

* 20 simulated networks ×
* 2 rooting settings (rooted/unrooted) ×
* 2 alignment types (full/third).

### Consensus Networks

The directory also contains the consensus networks produced using PhyloSummaries.jl.

There are four consensus files corresponding to the four analysis settings:

* rooted + full alignment,
* rooted + one-third alignment,
* unrooted + full alignment,
* unrooted + one-third alignment.

For example,

* `Consensus_nooutgroup_third.txt`

contains the consensus networks for the unrooted analyses performed on the one-third alignments.

Each consensus file contains 20 Newick networks, one per line. Line *i* corresponds to the consensus network for simulated network *i*.

### HMS and HWC Consensus Networks

For every simulated network and each analysis setting, the directory also contains the consensus networks obtained using the HMS and HWC methods.

Examples include:

* `HMS_outgroup_full_NetX.txt`
* `HWC_outgroup_full_NetX.txt`

where `X` denotes the simulated network.

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
* the path to the Julia executable,

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

Runs Squirrel on the full concatenated alignments while providing the root, producing both the inferred network and 100 bootstrap networks.

#### `4SQUIRREL_rooted_partial.py`

Runs Squirrel on the one-third alignments while providing the root.

#### `5SQUIRREL_unrooted_full.py`

Runs Squirrel on the full concatenated alignments without providing the root.

#### `6SQUIRREL_unrooted_partial.py`

Runs Squirrel on the one-third alignments without providing the root.

#### `7Huber_HWC_cmdlines.txt`

Contains the terminal commands required to compute the HMS and HWC consensus networks for all simulation settings.

#### `8PhyloSummaries.jl`

Julia implementation used to compute the consensus networks introduced in this work from the Squirrel bootstrap networks.

### Analysis and Figure Generation (Files 9–15)

The remaining files process the simulation results and reproduce the tables and figures presented in the paper.

* `9Process_edges_outgroup.R` – Processes edge-based metrics for the rooted analyses.
* `10Process_edges_nooutgroup.R` – Processes edge-based metrics for the unrooted analyses.
* `11Process_blobs_outgroup.R` – Processes blob, hybrid, circular order-based metrics for the rooted analyses.
* `12Process_blobs_nooutgroup.R` – Processes blob, hybrid, circular order-based metrics for the unrooted analyses.
* `13CreateExcels.R` – Compiles the processed results into summary spreadsheets.
* `14CreatePlots.R` – Generates the main figures presented in the paper.
* `15CreatePlots_disttoinput.R` – Generates figures summarizing the distance of the inferred consensus networks to the true simulated networks.


