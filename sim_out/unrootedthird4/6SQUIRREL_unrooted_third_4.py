"""
Bootstrapping consensus workflow for `physquirrel==1.0.6`.

Since `physquirrel==1.0.6` does not support
gene-block bootstrapping out of the box, this file provides a small wrapper
(`MSAWithGenes`) that adds:

- `gene_lengths` bookkeeping
- `bootstrap_genes(...)` sampling whole genes with replacement

Notes
-----
The workflow expects `.fasta` alignments and a text file of gene lengths.
When an outgroup is specified, `physquirrel` returns a directed network;
otherwise it returns a semi-directed network which is rooted at a leaf that
admits rooting.
"""

from __future__ import annotations

import os
import random
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, TypeVar


import physquirrel as psq
from pathlib import Path

T = TypeVar("T")

########################################################
# Newick string function (with internal node labels)
########################################################

def create_enewick_internal(net: psq.DirectedNetwork) -> str:
    """Create an eNewick string with stable internal node labels.

    This is a lightly adapted variant of `DirectedNetwork.create_enewick()` that
    labels internal nodes deterministically as `Node1`, `Node2`, ... (in the
    order returned by `net.internal_nodes()`).

    Parameters
    ----------
    net
        A directed phylogenetic network as returned by `physquirrel`.

    Returns
    -------
    str
        eNewick string (terminated with `;`).
    """
    internal_node_dict = {
        node: f"Node{str(i + 1)}" for i, node in enumerate(net.internal_nodes())
    }

    edge_list = list(net.edges)
    digraph: Dict[str, List[str]] = {}
    net._visited = {}  # type: ignore[attr-defined]
    in_degree: Dict[str, int] = {}

    for arc in edge_list:
        if arc[0] not in digraph:
            digraph[arc[0]] = []
        net._visited[arc[0]] = 0  # type: ignore[index]
        net._visited[arc[1]] = 0  # type: ignore[index]

        in_degree[arc[1]] = in_degree.get(arc[1], 0) + 1
        in_degree.setdefault(arc[0], 0)
        digraph[arc[0]].append(arc[1])

    root = "r"
    for node, degree in in_degree.items():
        if degree == 0:
            root = node

    net._reticulationVerticesFound = 0  # type: ignore[attr-defined]
    net._reticulationVertices = {}  # type: ignore[attr-defined]

    converted = _enewick_rec(
        net,
        source=root,
        digraph=digraph,
        in_degree=in_degree,
        internal_nodes_display=True,
        internal_node_dict=internal_node_dict,
    )
    return converted + ";"


def _enewick_rec(
    net: psq.DirectedNetwork,
    source: str,
    digraph: Mapping[str, Sequence[str]],
    in_degree: Mapping[str, int],
    internal_nodes_display: bool,
    internal_node_dict: Optional[Mapping[str, str]] = None,
) -> str:
    """Recursive helper for eNewick construction.

    Parameters
    ----------
    net
        Directed network being traversed.
    source
        Current node label.
    digraph
        Mapping from node -> list of children.
    in_degree
        In-degree map computed from `net.edges`.
    internal_nodes_display
        If True, internal nodes are labelled; if False only leaves are labelled.
    internal_node_dict
        Optional mapping to rename internal nodes.

    Returns
    -------
    str
        eNewick fragment for the subtree/network below `source`.
    """
    enewick = ""

    if in_degree[source] > 1:
        if source not in net._reticulationVertices:  # type: ignore[attr-defined]
            net._reticulationVerticesFound += 1  # type: ignore[attr-defined]
            ret_num = net._reticulationVerticesFound  # type: ignore[attr-defined]
            net._reticulationVertices[source] = ret_num  # type: ignore[attr-defined]
        else:
            ret_num = net._reticulationVertices[source]  # type: ignore[attr-defined]

    if net._visited[source] == 0:  # type: ignore[attr-defined]
        net._visited[source] = 1  # type: ignore[attr-defined]
        if source in digraph:
            enewick = "("
            for i, child in enumerate(digraph[source]):
                if i > 0:
                    enewick += ","
                enewick += _enewick_rec(
                    net,
                    source=child,
                    digraph=digraph,
                    in_degree=in_degree,
                    internal_nodes_display=internal_nodes_display,
                    internal_node_dict=internal_node_dict,
                )
            enewick += ")"

    if in_degree[source] > 1:
        enewick += f"#H{ret_num}:1.0"
    elif internal_nodes_display or source not in digraph:
        if internal_node_dict is not None and source not in net.leaves:
            enewick += f"{internal_node_dict[source]}:1.0"
        else:
            enewick += f"{source}:1.0"

    return enewick


########################################################
# Gene length reading function.
########################################################

def read_gene_lengths(gene_length_file: str) -> List[int]:
    import csv
    gene_lengths: List[int] = []
    with open(gene_length_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                gene_lengths.append(int(row[1]))  # second column is the length
            except ValueError:
                continue  # skips the header row ",x"
    return gene_lengths


########################################################
# Gene bootstrapping functionality
########################################################

def _bootstrap_sequences_by_genes(
    sequences: Mapping[str, str],
    gene_lengths: Sequence[int],
    *,
    seed: Optional[int] = None,
) -> Tuple[Dict[str, str], List[int]]:
    """Bootstrap an alignment by sampling whole genes with replacement.

    Parameters
    ----------
    sequences
        Mapping from taxon name to aligned sequence string.
    gene_lengths
        Lengths of consecutive gene blocks that partition each sequence.
        The sum must equal the sequence length.
    seed
        Optional random seed for reproducibility.

    Returns
    -------
    (dict[str, str], list[int])
        Bootstrapped sequences and the sampled gene-length list (one entry per
        sampled gene, with replacement).
    """
    if not sequences:
        return {}, []

    lengths = {len(s) for s in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("All sequences must have the same length.")
    (seq_len,) = tuple(lengths)
    if sum(gene_lengths) != seq_len:
        raise ValueError("Sum of gene lengths must equal sequence length.")

    rng = random.Random(seed)
    gene_indices = [rng.randrange(len(gene_lengths)) for _ in range(len(gene_lengths))]
    sampled_gene_lengths = [int(gene_lengths[i]) for i in gene_indices]

    # Precompute gene boundaries to avoid repeated summations.
    starts: List[int] = []
    pos = 0
    for gl in gene_lengths:
        starts.append(pos)
        pos += gl

    out: Dict[str, str] = {}
    for taxon, original_seq in sequences.items():
        parts: List[str] = []
        for gene_idx in gene_indices:
            start = starts[gene_idx]
            end = start + gene_lengths[gene_idx]
            parts.append(original_seq[start:end])
        out[taxon] = "".join(parts)
    return out, sampled_gene_lengths


class MSAWithGenes(psq.MSA):
    """`physquirrel.MSA` extended with `gene_lengths` and gene bootstrapping.

    Parameters
    ----------
    arg
        Same as `physquirrel.MSA`: a file path (`.fa/.fasta/.nex/.nexus`) or a
        dictionary mapping taxa to sequences.
    gene_lengths
        Optional list of gene lengths whose sum equals the sequence length.
        If provided, enables `bootstrap_genes`.
    """

    gene_lengths: Optional[List[int]]

    def __init__(self, arg: object, gene_lengths: Optional[Sequence[int]] = None):
        super().__init__(arg)

        lengths = {len(seq) for seq in self.sequences.values()}  # type: ignore[attr-defined]
        if len(lengths) != 1:
            if sum(gene_lengths_list) != seq_len:
                raise ValueError(
                  f"Sum of gene lengths ({sum(gene_lengths_list)}) must equal "
                  f"sequence length ({seq_len}). Gene lengths parsed: {gene_lengths_list[:10]}..."
                )
                raise ValueError("All sequences must have the same length.")
        (seq_len,) = tuple(lengths)

        if gene_lengths is None:
            self.gene_lengths = None
            return

        gene_lengths_list = list(map(int, gene_lengths))
        if sum(gene_lengths_list) != seq_len:
            
            raise ValueError(
                f"Sum of gene lengths ({sum(gene_lengths_list)}) must equal "
                f"sequence length ({seq_len}). Gene lengths parsed: {gene_lengths_list[:10]}..."
            )
            raise ValueError("Sum of gene lengths must equal sequence length.")
        self.gene_lengths = gene_lengths_list

    def bootstrap_genes(self, seed: Optional[int] = None) -> "MSAWithGenes":
        """Bootstrap by sampling whole gene blocks (requires `gene_lengths`).

        Parameters
        ----------
        seed
            Optional random seed for reproducibility.

        Returns
        -------
        MSAWithGenes
            A new MSA consisting of concatenated sampled gene blocks. The
            returned MSA keeps the same `gene_lengths`.
        """
        if self.gene_lengths is None:
            raise ValueError(
                "gene_lengths were not provided at initialization; cannot bootstrap by genes."
            )

        boot, sampled_gene_lengths = _bootstrap_sequences_by_genes(
            self.sequences,  # type: ignore[attr-defined]
            self.gene_lengths,
            seed=seed,
        )
        return MSAWithGenes(boot, gene_lengths=sampled_gene_lengths)


#############################################
# Processing a single alignment file.
#############################################



def process_file(
    args: Tuple[
        str,
        str,
        str,
        Sequence[int],
        int,
        Optional[str],
    ],
) -> None:
    """Process a single alignment file and write bootstrap networks.

    Parameters
    ----------
    args
        Tuple with entries

        - `file` : str
            Filename (must end with `.fasta` to be processed).
        - `folder` : str
            Input folder containing `file`.
        - `output_folder` : str
            Output folder where the per-file Newick output is appended.
        - `gene_lengths` : Sequence[int]
            Gene-length partition of the alignment. Must sum to the alignment
            length.
        - `nr_bootstraps` : int
            Number of bootstrap replicates. The unbootstrapped alignment is
            always included, so `nr_bootstraps + 1` networks are written.
        - `outgroup` : str | None
            If provided, networks are returned as directed networks rooted at
            this outgroup via `Q.squirrel(outgroup=...)`. If None, a rootable
            leaf is selected and the semi-directed network is converted to a
            directed one for eNewick output.

    Returns
    -------
    None
        This function is executed for its side effects: it appends one eNewick
        line per bootstrap replicate to `network_newick_<stem>.txt`.
    """
    file, folder, output_folder, gene_lengths, nr_bootstraps, outgroup = args

    if not file.endswith(".fasta"):
        return

    msa = MSAWithGenes(os.path.join(folder, file), gene_lengths=gene_lengths)

    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f"unrooted_network_newick_{file.replace('.fasta','')}.txt")

    for i in range(nr_bootstraps + 1):
        if i == 0:
            bootstrapped = msa
        else:
            bootstrapped = msa.bootstrap_genes(seed=i + 1)

        Q = bootstrapped.delta_heuristic(lam=0.15)

        if outgroup is None:
            networks, _ = Q.squirrel(include_score=True, triangles=False, all_networks=True)
        else:
            networks, _ = Q.squirrel(
                include_score=True,
                triangles=False,
                all_networks=True,
                outgroup=outgroup,
            )

        network = networks[0]

        if outgroup is None:
            newick1: Optional[str] = None
            for root in network.leaves:
                if network.is_rootable_at(root):
                    rooted_network = network.to_directed_network(root)
                    newick1 = create_enewick_internal(rooted_network)
                    break
            if newick1 is None:
                raise RuntimeError("Failed to find a valid rooting leaf.")
        else:
            assert isinstance(
                network, psq.DirectedNetwork
            ), "Expected a directed network when outgroup is specified"
            newick1 = create_enewick_internal(network)

        with open(output_file, "a", encoding="utf-8") as f:
            f.write(newick1 + "\n")

    print(f"Completed file: {file}")


########################################################
# Main function (includes parallel processing)
########################################################

def main(
    nr_bootstraps: int,
    outgroup: Optional[str],
    fasta_folder: str,
    output_folder: str,
    gene_length_file: str,
) -> None:
    """Run the bootstrapping workflow with multiprocessing.

    Parameters
    ----------
    nr_bootstraps
        Number of bootstrap replicates (in addition to the original MSA, which
        is always included).
    outgroup
        If provided, root at this outgroup via `Q.squirrel(outgroup=...)`.
    fasta_folder
        Folder containing `.fasta` files.
    output_folder
        Folder to write `unrooted_network_newick_*.txt` files into.
    gene_length_file
        Path to a gene-length text file.
    """
    gene_lengths = read_gene_lengths(gene_length_file)
    sorted_files = sorted(os.listdir(fasta_folder))

    args_list = [
        (file, fasta_folder, output_folder, gene_lengths, nr_bootstraps, outgroup)
        for file in sorted_files
        if file.endswith(".fasta")
    ]

    with Pool(processes=cpu_count()) as pool:
        pool.map(process_file, args_list)




# Root directory (change only this)
base_dir = Path(
    "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-third10bf"
).expanduser()

# Derived paths
sim_dir = base_dir / "sim_out/unrootedthird4"
fasta_folder = sim_dir
gene_length_file = sim_dir / "seq_lens_third.csv"
output_folder = sim_dir

# Other parameters
nr_bootstraps = 100
outgroup = None  # or None

if __name__ == "__main__":
    main(
        nr_bootstraps=nr_bootstraps,
        outgroup=outgroup,
        fasta_folder=str(fasta_folder),
        output_folder=str(output_folder),
        gene_length_file=str(gene_length_file),
    )
