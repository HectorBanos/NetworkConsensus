#!/usr/bin/env bash
#
# ML trees with bootstrap support for the two HIV alignments in data/.
#
# Settings: GTR+GAMMA, ML search + bootstrapping in one run (--all), autoMRE
# bootstopping (max 1000 replicates), seed 42, --threads auto. Support values
# are Felsenstein bootstrap proportions (FBP) on the best-scoring ML tree.
#
# A second step derives the extended majority-rule (MRE) consensus of the
# bootstrap replicate trees (--consense MRE); its node labels are bootstrap
# split frequencies, not FBP, and it carries no branch lengths.
#
# Version used: raxml-ng v2.0.2 (conda/miniforge3), on Debian 6.1.
#
# Per-run raxml-ng files go to trees_raxmlng/logs/<name>/; each summary block
# holds four trees (ML best and MRE consensus, with and without support) plus
# the run details, collected into trees_raxmlng/raxml_trees.txt.

set -euo pipefail

RAXML=~/miniforge3/bin/raxml-ng

# Paths are relative to this script, so the folder can be moved or cloned.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$SCRIPT_DIR/data"
OUTPUT_DIR="$SCRIPT_DIR/trees_raxmlng"
SUMMARY="$OUTPUT_DIR/raxml_trees.txt"
mkdir -p "$OUTPUT_DIR"

# raxml-ng writes ~9 files per run; keep them, one subfolder per alignment.
LOGS_DIR="$OUTPUT_DIR/logs"

ALIGNMENTS=(references_all references_noJH)

: > "$SUMMARY"

for name in "${ALIGNMENTS[@]}"; do
    echo "=== $name ==="
    mkdir -p "$LOGS_DIR/$name"
    prefix="$LOGS_DIR/$name/$name"

    "$RAXML" --all \
        --msa "$INPUT_DIR/$name.fasta" \
        --model GTR+G \
        --prefix "$prefix" \
        --bs-trees autoMRE \
        --seed 42 \
        --threads auto \
        --redo

    # Extended majority-rule consensus of the bootstrap replicate trees. Written
    # to a temp prefix so only the ML run's files remain under logs/<name>/; the
    # tree itself is collected into the summary below.
    mre_tmp="$(mktemp -d)"
    "$RAXML" --consense MRE \
        --tree "$prefix.raxml.bootstraps" \
        --prefix "$mre_tmp/mre" \
        --redo
    mre_tree="$mre_tmp/mre.raxml.consensusTreeMRE"

    # Collect the results: ML tree with FBP support and the MRE bootstrap
    # consensus (each with and without support), plus the run details that
    # make the numbers interpretable (taxa/sites, likelihood, bootstraps used).
    {
        echo "################################################################"
        echo "# $name"
        echo "################################################################"
        echo
        echo "Alignment: input/$name.fasta"
        grep -E "^\[.*\] Loaded alignment|Alignment comprises|^Gaps:|^Invariant sites:" \
            "$prefix.raxml.log" || true
        grep -E "Fully undetermined columns" "$prefix.raxml.log" || true
        echo
        echo "Final ML log-likelihood: $(grep -oP 'Final LogLikelihood: \K.*' "$prefix.raxml.log" || echo 'n/a')"
        echo "Bootstrap replicates: $(grep -oP 'Bootstrapping converged after \K[0-9]+' "$prefix.raxml.log" || echo 'did not converge (max used)')"
        echo "Model: $(grep -oP '^   Substitution model: \K.*' "$prefix.raxml.log" | head -1 || echo 'GTR+G')"
        echo
        echo "--- Best ML tree (with FBP support) ---"
        cat "$prefix.raxml.support"
        echo
        echo "--- Best ML tree (no support values) ---"
        cat "$prefix.raxml.bestTree"
        echo
        echo "--- MRE consensus tree (with support) ---"
        cat "$mre_tree"
        echo
        echo "--- MRE consensus tree (no support values) ---"
        # Strip the split-frequency labels that follow each closing paren.
        sed -E 's/\)[0-9]+/)/g' "$mre_tree"
        echo
    } >> "$SUMMARY"

    rm -rf "$mre_tmp"
done

echo
echo "Wrote $SUMMARY"
