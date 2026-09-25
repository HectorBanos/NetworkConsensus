"""Infer a phylogenetic network per sampled taxon combination.

For each replicate we pick one representative sequence per subtype, cut the
corresponding sub-MSA out of the full alignment, and run squirrel on it. We do
this in two rooting modes: rooted at an outgroup (a DirectedNetwork) and
unrooted (a SemiDirectedNetwork, no outgroup). Each inferred network is written
twice: once with the full FASTA taxa names, and once with those names replaced
by their subtype labels.

Pipeline: sample taxon combinations -> extract sub-MSA -> infer network ->
convert taxa names to subtypes -> write both eNewick variants.

Layout: the input alignment lives in 'data/', results are written to
'networks_squirrel/'.

Requires (version used): Python 3.11, physquirrel 1.0.9
(pip install physquirrel==1.0.9).
"""

import os
import re
import random
from math import prod
from importlib.metadata import version
from typing import List, Dict

import physquirrel as psq


#### Subtype data ############################################################
# Subtype-wise HIV sequence dictionary derived from RIP FASTA headers. Each
# subtype maps to the sequences that represent it; one is sampled per replicate.

RIP_SUBTYPE_SEQUENCES = {
    'A1': [
        'A1.KE.94.Q23_17.AF004885',
        'A1.SE.94.SE7253.AF069670',
        'A1.UG.92.92UG037.U51190',
        'A1.UG.98.98UG57136.AF484509',
    ],
    'A2': [
        'A2.CD.97.97CDKTB48.AF286238',
        'A2.CM.01.01CM_1445MV.GU201516',
        'A2.CY.94.94CY017_41.AF286237',
    ],
    'A3': [
        'A3.SN.01.DDI579.AY521629',
        'A3.SN.01.DDJ369.AY521631',
        'A3.SN.96.DDJ360.AY521630',
    ],
    'A4': [
        'A4.CD.97.97CD_KCC2.AM000053',
        'A4.CD.97.97CD_KTB13.AM000054',
    ],
    'A6': [
        'A6.BY.13.PV85.KT983615',
        'A6.GE.99.99GEMZ011.DQ207944',
        'A6.IT.02.60000.EU861977',
        'A6.UA.00.98UA0116.AF413987',
    ],
    'B': [
        'B.FR.83.HXB2_LAI_IIIB_BRU.K03455',
        'B.NL.99.671_99T12.AY423381',
        'B.TH.90.BK132.AY173951',
        'B.US.98.1058_11.AY331295',
    ],
    'C': [
        'C.BR.92.BR025_d.U52953',
        'C.ET.86.ETH2220.U46016',
        'C.IN.95.95IN21068.AF067155',
        'C.ZA.04.04ZASK146.AY772699',
    ],
    'D': [
        'D.CD.83.ELI.K03454',
        'D.CM.01.01CM_4412HAL.AY371157',
        'D.TZ.01.A280.AY253311',
        'D.UG.94.94UG114.U88824',
    ],
    'F1': [
        'F1.BE.93.VI850.AF077336',
        'F1.BR.93.93BR020_1.AF005494',
        'F1.FI.93.FIN9363.AF075703',
        'F1.FR.96.96FR_MP411.AJ249238',
    ],
    'F2': [
        'F2.CM.02.02CM_0016BBY.AY371158',
        'F2.CM.95.95CM_MP255.AJ249236',
        'F2.CM.95.95CM_MP257.AJ249237',
        'F2.CM.97.CM53657.AF377956',
    ],
    'G': [
        'G.BE.96.DRCBL.AF084936',
        'G.KE.93.HH8793_12_1.AF061641',
        'G.NG.92.92NG083_JV10832.U88826',
        'G.SE.93.SE6165_G6165.AF061642',
    ],
    'H': [
        'H.BE.93.VI991.AF190127',
        'H.BE.93.VI997.AF190128',
        'H.CF.90.056.AF005496',
        'H.GB.00.00GBAC4001.FJ711703',
    ],
    'J': [
        'J.CM.04.04CMU11421.GU237072',
        'J.SE.93.SE9280_7887.AF082394',
        'J.SE.94.SE9173_7022.AF082395',
    ],
    'K': [
        'K.CD.97.97ZR_EQTB11.AJ249235',
        'K.CM.96.96CM_MP535.AJ249239',
    ],
    'N': [
        'N.CM.02.DJO0131.AY532635',
        'N.CM.04.04CM_1015_04.DQ017382',
        'N.CM.95.YBF30.AJ006022',
        'N.CM.97.YBF106.AJ271370',
    ],
    # Single-element group: the hybrid taxon, always present in every sample.
    'A6B': [
        '03_A6B.RU.97.KAL153_2.AF193276',
    ],
}


#### Sampling ################################################################
# Each subtype contributes exactly one sequence per replicate. A replicate is
# therefore one element of the Cartesian product of the per-subtype groups.
# Rather than build that (astronomically large) product, we sample integers in
# [0, total) and decode each one into a combination.

def sample_uniform_combinations(sequence_groups, sample_size, seed=0):
    """Uniformly sample unique tuples from a Cartesian product of groups.

    Returns one tuple per replicate, holding one sequence per group in group
    order. `sample_size` is capped at the number of combinations that exist.
    """
    rng = random.Random(seed)

    # Total size of the Cartesian product. rng.sample draws without replacement,
    # so every returned combination is unique.
    group_sizes = [len(group) for group in sequence_groups]
    total_combinations = prod(group_sizes)
    k = min(sample_size, total_combinations)
    sampled_indices = rng.sample(range(total_combinations), k=k)

    # Decode each sampled integer into one concrete combination.
    sampled_combinations = []
    for idx in sampled_indices:
        combination = []
        for size, group in zip(reversed(group_sizes), reversed(sequence_groups)):
            idx, choice_idx = divmod(idx, size)
            combination.append(group[choice_idx])
        # Built back-to-front above, so reverse to restore group order.
        sampled_combinations.append(tuple(reversed(combination)))

    return sampled_combinations


#### Sub-MSA extraction ######################################################
# squirrel runs on one small alignment per replicate, not the full FASTA, so we
# slice out just the sampled taxa.

def extract_subMSA(msa: psq.MSA, taxa_list: List[str]) -> psq.MSA:
    """Return a new psq.MSA containing only the taxa in `taxa_list`, in that order.

    Raises ValueError if any requested taxon is missing from `msa`.
    """
    if not hasattr(msa, "sequences") or not hasattr(msa, "taxa"):
        raise TypeError("`msa` must be a physquirrel.MSA instance (or have `.sequences` and `.taxa`).")
    missing = [t for t in taxa_list if t not in msa.taxa]
    if missing:
        raise ValueError(f"The following taxa are not present in the input MSA: {missing}")

    # Preserve the provided order in taxa_list
    sequences: Dict[str, str] = {t: msa.sequences[t] for t in taxa_list}

    return psq.MSA(sequences)


#### Newick conversion #######################################################
# Inferred networks carry full FASTA names (e.g. 'B.TH.90.BK132.AY173951').
# For output for which we can run a consensus we also emit a copy labelled by subtype (e.g. 'B').

def build_taxa_to_subtype_map(subtypes=None):
    """Map full taxa names to their subtypes, for the given subtypes (default: all)."""
    if subtypes is None:
        subtypes = list(RIP_SUBTYPE_SEQUENCES.keys())

    taxa_to_subtype = {}
    for subtype in subtypes:
        for taxa in RIP_SUBTYPE_SEQUENCES[subtype]:
            taxa_to_subtype[taxa] = subtype

    return taxa_to_subtype


def replace_taxa_in_newick(newick_string, taxa_to_subtype):
    """Replace all taxa names in a Newick string with their subtypes."""
    result = newick_string

    # Sort by length descending to replace longer names first (avoid partial replacements)
    sorted_taxa = sorted(taxa_to_subtype.keys(), key=len, reverse=True)

    for taxa in sorted_taxa:
        subtype = taxa_to_subtype[taxa]
        # Use word boundaries to avoid partial replacements
        result = re.sub(r'\b' + re.escape(taxa) + r'\b', subtype, result)

    return result


#### Pipeline ################################################################

def sample_squirrel_networks(
    fasta_file,
    output_folder,
    output_name,
    outgroup=None,
    subtypes=None,
    n_replicates=100,
    lam=0.3,
    seed=0,
):
    """Sample taxon combinations, infer a network for each, and write both eNewick variants.

    Args:
        fasta_file: path to the input FASTA alignment.
        output_folder: directory for the results. Created if it does not exist.
        output_name: base name for the outputs. Two files are written:
            '<output_name>_full.txt' (full taxa names) and
            '<output_name>_subtypes.txt' (subtype labels).
        outgroup: subtype used to root each network (e.g. 'N'), or None for the
            unrooted/semi-directed case. If given it must be in `subtypes`; if
            None, squirrel returns a SemiDirectedNetwork and create_enewick()
            serialises it with an arbitrary rooting.
        subtypes: which subtypes from RIP_SUBTYPE_SEQUENCES to include, one
            sequence sampled per subtype per replicate. Defaults to all of them.
        n_replicates: how many taxon combinations to sample.
        lam: lambda for the delta heuristic.
        seed: seed for reproducible sampling.

    Returns:
        A (raw_path, converted_path) tuple of the two files written.
    """
    # Record which physquirrel produced these networks: inference behaviour can
    # change between versions, so the version is part of the result.
    print(f'Using physquirrel version: {version("physquirrel")}')

    # Default to every subtype in the dictionary, in declaration order.
    if subtypes is None:
        subtypes = list(RIP_SUBTYPE_SEQUENCES.keys())

    # Validate up front.
    unknown = [s for s in subtypes if s not in RIP_SUBTYPE_SEQUENCES]
    if unknown:
        raise ValueError(f"Unknown subtypes (not in RIP_SUBTYPE_SEQUENCES): {unknown}")
    if outgroup is not None and outgroup not in subtypes:
        raise ValueError(f"Outgroup {outgroup!r} is not among the selected subtypes: {subtypes}")

    # Combinations follow the order of `subtypes`, so the outgroup's subtype
    # position tells us which taxon in each combination is the outgroup. In the
    # unrooted case there is no outgroup and hence no such index.
    sequence_groups = [RIP_SUBTYPE_SEQUENCES[subtype] for subtype in subtypes]
    outgroup_index = subtypes.index(outgroup) if outgroup is not None else None

    msa = psq.MSA(fasta_file)
    sampled_combinations = sample_uniform_combinations(
        sequence_groups, sample_size=n_replicates, seed=seed
    )
    taxa_to_subtype = build_taxa_to_subtype_map(subtypes)

    # Two parallel outputs: same networks, different leaf labels. The seqname
    # (full FASTA names) versions go in a 'seqname/' subfolder; the subtypename
    # versions stay in output_folder, since that is what the consensus methods read.
    seqname_folder = os.path.join(output_folder, 'seqname')
    os.makedirs(seqname_folder, exist_ok=True)
    raw_path = os.path.join(seqname_folder, f'{output_name}_seqname.txt')
    converted_path = os.path.join(output_folder, f'{output_name}_subtypename.txt')

    total = len(sampled_combinations)

    with open(raw_path, 'w', encoding='utf-8') as raw_handle, \
         open(converted_path, 'w', encoding='utf-8') as converted_handle:

        for count, combination in enumerate(sampled_combinations, start=1):
            # Infer the best network for this replicate: slice the alignment,
            # score it with the delta heuristic, then let squirrel build the
            # network. With an outgroup squirrel returns a network rooted at this
            # replicate's outgroup sequence; without one it returns a
            # semi-directed network.
            sub_msa = extract_subMSA(msa, list(combination))
            Q = sub_msa.delta_heuristic(lam=lam, weight=True)
            if outgroup is None:
                best_network = Q.squirrel()
            else:
                outgroup_taxon = combination[outgroup_index]
                best_network = Q.squirrel(outgroup=outgroup_taxon)

            # Write the network twice, keeping the two files line-for-line aligned.
            # create_enewick() on a semi-directed network uses an arbitrary rooting.
            enewick = best_network.create_enewick()
            raw_handle.write(enewick + '\n')
            converted_handle.write(replace_taxa_in_newick(enewick, taxa_to_subtype) + '\n')

            # Overwrite the same terminal line rather than scrolling.
            print(f'\rNetworks inferred: {count}/{total}', end='', flush=True)

    print(f'\rNetworks inferred: {total}/{total}')
    print(f'Wrote {raw_path}')
    print(f'Wrote {converted_path}')

    return raw_path, converted_path


#### Configuration ###########################################################
# Shared settings for every run. Paths are relative to this script, so the
# folder can be moved or cloned without editing them.

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FASTA_FILE = os.path.join(SCRIPT_DIR, 'data', 'HIV1_RIP_2023_genome_DNA.fasta')
OUTPUT_FOLDER = os.path.join(SCRIPT_DIR, 'networks_squirrel')
OUTGROUP = 'N'
N_REPLICATES = 1000
LAM = 0.3
SEED = 0

# Two subtype selections: everything, and everything except J and H.
ALL_SUBTYPES = ['A1', 'A2', 'A3', 'A4', 'A6', 'B', 'C', 'D', 'F1', 'F2',
                'G', 'H', 'J', 'K', 'N', 'A6B']
WITHOUT_JH = ['A1', 'A2', 'A3', 'A4', 'A6', 'B', 'C', 'D', 'F1', 'F2',
              'G', 'K', 'N', 'A6B']


def run_mode(rooting):
    """Sample networks for both subtype selections in one rooting mode.

    `rooting` is 'rooted' (each network rooted at the outgroup subtype OUTGROUP)
    or 'unrooted' (semi-directed, no outgroup). It also names the outputs, e.g.
    'squirrel_rooted_sampling_HIV_all' / 'squirrel_unrooted_sampling_HIV_noJH'.
    """
    if rooting not in ('rooted', 'unrooted'):
        raise ValueError("rooting must be 'rooted' or 'unrooted'")
    outgroup = OUTGROUP if rooting == 'rooted' else None

    for selection, subtypes in (('all', ALL_SUBTYPES), ('noJH', WITHOUT_JH)):
        sample_squirrel_networks(
            fasta_file=FASTA_FILE,
            output_folder=OUTPUT_FOLDER,
            output_name=f'squirrel_{rooting}_sampling_HIV_{selection}',
            outgroup=outgroup,
            subtypes=subtypes,
            n_replicates=N_REPLICATES,
            lam=LAM,
            seed=SEED,
        )


if __name__ == '__main__':
    # Full analysis: both rooting modes x both subtype selections -> 8 files.
    run_mode('rooted')
    run_mode('unrooted')
