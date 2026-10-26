#!/usr/bin/perl 

# HWC-consensus.pl
# Date: 2026-05-15
#
# ------------------------------------------------------------------------
# Description:
#   Computes the hardwired cluster (HWC) consensus phylogenetic network
#   as described in Nielsen et al. (2023), given a set of input networks.
#
# Input:
#   A CSV-like text file containing extended Newick strings.
#   Each entry may optionally include a weight.
#
# Output:
#   A single consensus network represented as an extended Newick string.
#
# Consensus Modes:
#   --mode strict     : Includes only clusters present in all input networks
#   --mode majority   : Includes clusters present in >50% of networks
#   --mode custom     : Includes clusters meeting a user-defined threshold
#
# Parameters:
#   --t <float>       : Threshold in (0,1), required for custom mode.
#                      A cluster is included if it appears in at least
#                      this fraction of the input networks.
#
# Example Usage:
#   perl HWC-consensus.pl \
#       --in=input-networks.txt \
#       --out=consensus-net.txt \
#       --mode=custom \
#       --t=0.75
#
# Notes:
#   - Input networks must be valid extended Newick format.
#   - Weights (if provided) should be numeric and associated per network.
#   - The script assumes consistent taxon labels across all input networks.
# ------------------------------------------------------------------------

use strict;
use warnings;

# Command-line option parsing
use Getopt::Long;

# External modules:
# Bio::PhyloNetwork:
#   Provides parsing and manipulation of phylogenetic networks
use Bio::PhyloNetwork;

# Graph:
#   General-purpose graph structure used for representing clusters
#   and relationships during consensus construction
use Graph;


#######################################################################
# ------------------- COMMAND LINE HANDLING ---------------------------
#######################################################################

# Variables for command-line arguments
my ($input, $output, $mode, $threshold);

# Parse command-line options
# --in        : input file containing networks
# --out       : output file for consensus network (optional)
# --mode      : consensus mode (strict, majority, custom)
# --threshold : threshold for custom mode (float in [0,1])
GetOptions(
    'in=s'        => \$input,
    'out=s'       => \$output,
    'mode=s'      => \$mode,
    'threshold=f' => \$threshold
) or die "Error in command line arguments\n";

# Set defaults if not provided
$mode = 'majority' unless defined $mode;
$threshold = 0.5 unless defined $threshold;

# Validate threshold range
# (Used only in custom mode, but validated globally for simplicity)
die "Threshold must be between 0 and 1\n"
	unless $threshold >= 0 && $threshold <= 1;

# Ensure input file is provided
die "Missing network input filename.\n"
	unless defined $input;

# Attempt to open input file early to fail fast if it is invalid
open(my $in, "<", $input)
	or die "Could not open file '$input'. $!";

# Output handling:
# If no output file is specified, results will be printed to STDOUT
my $terminal = 0;

if (!defined $output){
	$terminal = 1;
}
# If output file exists, ask user how to proceed
elsif (-e $output){

	print "\nThe file $output already exists.\n",
			"Type OVERWRITE or APPEND:\n";

    # Read user response from STDIN
	my $answer = <STDIN>;
	chomp($answer);
	$answer = uc($answer);  # Normalize input for comparison

	if ($answer eq "OVERWRITE"){
		# Delete existing file before writing new output
		unlink($output)
		or die "Could not delete $output: $!";
	}
	elsif ($answer ne "APPEND"){
		# Reject any response other than OVERWRITE or APPEND
		die "Invalid response\n";
	}
}



#######################################################################
# ------------------- INPUT PARSING -----------------------------------
#######################################################################

# Read input file containing phylogenetic networks.
#
# Expected format (per line):
#   <extended Newick string> [, weight]
#
# Notes:
#   - The Newick string must end with a semicolon.
#   - The weight is optional; if omitted, a default weight of 1 is used.
#   - Lines that are empty or contain only whitespace are ignored.
#   - The first line is skipped (e.g., may represent a reference/"true" network).

open(my $fh, '<', $input)
    or die "Cannot open $input: $!\n";

my @networks;  # Array of hashrefs: { newick => ..., weight => ... }

# Discard first line (if it is not part of the network sample, for example, if it is the target network)
# <$fh>;

while (my $line = <$fh>) {
	chomp $line;

	# Trim leading and trailing whitespace
	$line =~ s/^\s+|\s+$//g;

	# Skip empty lines
	next if $line eq '';

	my ($newick, $weight);

	# Parse line into Newick string and optional weight
    if ($line =~ /^(.*?;)\s*(?:,\s*([0-9.eE+-]+)\s*)?$/) {
		$newick = $1;

		# Assign default weight = 1 if not provided
		$weight = defined $2 ? $2 : 1;
	}
	else {
		# Fail early if line does not match expected format
		die "Line does not look like extended Newick:\n$line\n";
	}

	# Store parsed network as a hash reference
	push @networks, {
		newick => $newick,
		weight => $weight,
	};

	# Debug (optional): track number of parsed networks
	# print "Number of networks: " . scalar(@networks) . "\n";
}

close($fh);



#######################################################################
# -------- CONSENSUS HARDWIRED CLUSTERS COMPUTATION -------------------
#######################################################################

# Hash tracking total weighted frequency of each cluster across networks.
# Key   : canonical cluster representation (sorted taxa joined by comma)
# Value : sum of weights of networks in which the cluster appears
my %cluster_counts;

foreach my $net (@networks) {

	# Use provided weight or default to 1 if missing
	my $weight = $net->{weight} // 1;

	# Extract all clusters from the current network
	my @clusters = get_Clusters($net);

	# Debug output: list clusters found in current network
	# print "Clusters in current network:\n";
	# foreach my $cluster (@clusters) {
	#    print join(",", @$cluster), "\n";
	# }
	#print "\n\n";

	# Track clusters already counted within this network
	# Prevents inflating counts due to duplicate clusters in same network
	my %seen_keys;

	foreach my $cluster (@clusters) {

		# Convert cluster into canonical representation:
		# - sort taxa so ordering does not affect identity
		# - join into string key for hashing
		my @sorted = sort @$cluster;
		my $key = join(",", @sorted);

		# Skip duplicate occurrences within the same network
		next if $seen_keys{$key};

		# Accumulate weighted presence across networks
		$cluster_counts{$key} += $weight;

		# Mark cluster as seen in this network
		$seen_keys{$key} = 1;
	}
}


# Compute total weight across all networks.
# This is used later to normalize cluster frequencies (e.g., for thresholds)
my $total_weight = 0;
$total_weight += $_->{weight} // 1 for @networks;



#######################################################################
# ------------------- CONSENSUS FILTERING -----------------------------
#######################################################################

# Select clusters that satisfy the chosen consensus rule.
#
# cluster_counts:
#   weighted number of networks in which each cluster appears
#
# total_weight:
#   sum of all network weights (used for normalization)
#
# Result:
#   @consensus_clusters contains cluster keys that pass the filter

my @consensus_clusters;

foreach my $key (keys %cluster_counts) {

	# Convert weighted counts into relative frequency
	# (i.e., proportion of total network weight supporting this cluster)
	my $freq = $cluster_counts{$key} / $total_weight;

    # Apply consensus rule depending on selected mode
		if (keep_cluster($freq, $mode, $threshold)) {
			push @consensus_clusters, $key;
		}
}


# Debug output: final selected clusters
# print "Clusters in consensus network\n";
# foreach my $cluster (@consensus_clusters) {
#	print $cluster, "\n";
# }
# print "\n";



# ----------------------------------------------------------------------
# Decide whether a cluster should be included in the consensus network.
#
# Inputs:
#   $freq      : normalized support (0..1)
#   $mode      : consensus mode (strict | majority | custom)
#   $threshold : user-defined cutoff for custom mode
#
# Returns:
#   1 if cluster is included, 0 otherwise
# ----------------------------------------------------------------------

sub keep_cluster {
    my ($freq, $mode, $threshold) = @_;

    if ($mode eq 'strict') {

        # Strict consensus: cluster must appear in ALL networks
        return $freq == 1;

    }
    elsif ($mode eq 'majority') {

        # Majority rule: cluster must appear in more than half
        return $freq > 0.5;

    }
    else {  # custom threshold mode

        # User-defined cutoff (e.g., 0.7 means 70% support required)
        return $freq >= $threshold;
    }
}



#######################################################################
# -------- GRAPH (HASSE DIAGRAM) CONSTRUCTION -------------------------
#######################################################################

# Convert consensus cluster representation back into list-of-taxa format.
#
# Earlier step stored clusters as:
#   "A,B,C"
#
# For graph construction we need:
#   [ "A", "B", "C" ]
#
# This restores structured data for lattice/Hasse operations.

my @fixed_clusters;

foreach my $c (@consensus_clusters) {

    # Split canonical cluster string into individual taxa labels
    my @arr = split /,/, $c;

    push @fixed_clusters, \@arr;
}

# Build Hasse diagram (partial order graph) from cluster inclusion relations.
#
# The Hasse graph represents subset/superset relationships between clusters:
#   A ⊂ B edges define hierarchy without transitive edges
#
# Input:
#   reference to array of cluster sets (arrayrefs of taxa)
#
# Output:
#   graph object encoding partial order structure

my $g_hasse = build_hasse_graph(\@fixed_clusters);

# Post-processing step to uncompress the network (meaning that reticulations have out-degree one)
uncompress($g_hasse);



#######################################################################
# ------------------- NODE RELABELING ---------------------------------
#######################################################################

# Create a working copy of the Hasse graph so that node labels can be
# safely modified without affecting the original structure.
my $g = $g_hasse->deep_copy_graph;

# Counters for generating new node labels
my $i_count = 1;  # internal tree-like nodes
my $h_count = 1;  # reticulation (hybrid) nodes

# Mapping from original vertex labels to new labels
my %rename;

foreach my $v ($g->vertices) {

	my $out = $g->out_degree($v);
	my $in  = $g->in_degree($v);

	# Skip leaves (terminal taxa): these retain original labels
	next if $out == 0;

	# Reticulation nodes:
	# Nodes with multiple parents (in-degree > 1)
	# These represent hybridization / network events
	if ($in > 1) {
		$rename{$v} = "#H$h_count";
		$h_count++;
	}

	# Internal tree nodes:
	# Standard branching points in the hierarchy
	else {
		$rename{$v} = "I$i_count";
		$i_count++;
	}
}

# ----------------------------------------------------------------------
# Rebuild graph with renamed internal structure
# ----------------------------------------------------------------------

my $g2 = Graph->new(directed => 1);

# Add all vertices with updated labels (or original for leaves)
foreach my $v ($g->vertices) {
	my $new = $rename{$v} // $v;  # leaves remain unchanged
	$g2->add_vertex($new);
}

# Recreate edges using renamed endpoints
foreach my $e ($g->edges) {
	my ($src, $dst) = @$e;

	my $new_src = $rename{$src} // $src;
	my $new_dst = $rename{$dst} // $dst;

	$g2->add_edge($new_src, $new_dst);
}



#######################################################################
# ------------------- OUTPUT ------------------------------------------
#######################################################################

# Convert the final relabeled graph into an extended Newick representation.
#
# The Graph object ($g2) encodes the consensus network structure.
# Bio::PhyloNetwork is used to translate this graph into an
# extended Newick string suitable for downstream phylogenetic tools.

my $consensusnet = new Bio::PhyloNetwork(-graph => $g2);

# Generate full extended Newick representation of the network
# (includes internal nodes and reticulation structure)
my $str = $consensusnet->eNewick_full();


# Print to STDOUT or file depending on user input
my $OUT;

if ($terminal) {
	# No output file provided → write to STDOUT
	open($OUT, '>&', \*STDOUT)
	or die "Cannot open STDOUT: $!";
}
else {
	# Append results to output file
	open($OUT, '>>', $output)
	or die "Could not open file '$output': $!";
}

print $OUT $str . "\n";
close($OUT);



#######################################################################
# ------------------- CORE FUNCTION: GET HARDWIRED CLUSTERS ----------
#######################################################################

# Extract hardwired clusters from a phylogenetic network.
#
# A hardwired cluster for node u is defined as the set of all leaf
# descendants reachable from u in the directed network graph.
#
# Input:
#   $net : hash reference containing an extended Newick string
#
# Output:
#   list of clusters (each cluster is an array reference of leaf labels)

sub get_Clusters {

	my ($net) = @_;
	my $newick = $net->{newick};

	die "Newick string undefined!" unless defined $newick;

	# Parse extended Newick into a network object
	my $network = Bio::PhyloNetwork->new(-eNewick => $newick);
	my $graph = $network->graph or die "Graph undefined";

	# Identify leaf nodes (nodes with no outgoing edges)
	my @nodes = $graph->vertices;
	my %leaves = map { $_ => 1 }
		grep { $graph->out_degree($_) == 0 } @nodes;

	my @clusters;

	# For each node, compute its descendant leaf set
	foreach my $u (@nodes) {

		my %cluster;

		# Include node itself if it is a leaf
		$cluster{$u} = 1 if exists $leaves{$u};

		# Add all descendant leaves reachable from this node
		foreach my $v ($graph->all_successors($u)) {
			$cluster{$v} = 1 if exists $leaves{$v};
		}

		# Convert to sorted canonical representation
		my @leaf_set = sort keys %cluster;

		push @clusters, \@leaf_set;
	}

	return @clusters;
}



#######################################################################
# ------------------- CORE FUNCTION: HASSE DIAGRAM --------------------
#######################################################################

# Construct a Hasse diagram representing subset relations between clusters.
#
# Nodes represent clusters.
# Directed edges represent *cover relations* in the subset partial order:
#   A → B if A ⊂ B and no intermediate cluster C exists such that
#   A ⊂ C ⊂ B.
#
# Input:
#   reference to array of clusters (each cluster is arrayref of taxa)
#
# Output:
#   directed Graph object encoding the Hasse diagram

sub build_hasse_graph {

	my ($clusters_ref) = @_;
	my @clusters = @$clusters_ref;

	my $g = Graph->new(directed => 1);

	# Maps canonical cluster string → set representation
	my %set;
	my @keys;

	foreach my $c (@clusters) {

		my @sorted = sort @$c;
		my $key = join(",", @sorted);

		push @keys, $key;
		$set{$key} = { map { $_ => 1 } @sorted };

		$g->add_vertex($key);
	}

	# Sort clusters by size (bottom-up construction of partial order)
	@keys = sort {
		scalar(keys %{$set{$a}}) <=> scalar(keys %{$set{$b}})
	} @keys;

    # Check subset relation between two clusters
	sub is_subset {
		my ($A, $B) = @_;

		foreach my $x (keys %$A) {
			return 0 unless $B->{$x};
		}

		return 1;
	}

	# Build Hasse edges (transitive reduction of subset lattice)
	for (my $i = 0; $i < @keys; $i++) {

		my $A = $set{$keys[$i]};

		for (my $j = $i + 1; $j < @keys; $j++) {

			my $B = $set{$keys[$j]};

			# Skip if not a subset relation
			next unless is_subset($A, $B);

			my $covered = 1;

			# Check for intermediate clusters (transitive reduction)
			for (my $k = $i + 1; $k < $j; $k++) {

				my $C = $set{$keys[$k]};

				if (is_subset($A, $C) && is_subset($C, $B)) {
					$covered = 0;
					last;
				}
			}

			# Add edge only if it is a direct cover relation
			if ($covered) {
				$g->add_edge($keys[$j], $keys[$i]);
			}
		}
	}

	return $g;
}



#######################################################################
# ------------------- GRAPH CLEANUP -----------------------------------
#######################################################################

# Insert subdivision ("cut") nodes to expand compressed graph structure.
#
# Two cases handled:
#   1. Reticulation-like structure (indegree ≥2, outdegree ≥2)
#   2. Root/terminal compression case (indegree ≥2, outdegree = 0)

sub uncompress {

	my ($g) = @_;

	my @vertices = $g->vertices;

	foreach my $v (@vertices) {

		my @pred = $g->predecessors($v);
		my @succ = $g->successors($v);

		my $indegree  = scalar @pred;
		my $outdegree = scalar @succ;

		# Only modify structurally complex nodes
		next unless (
			($indegree >= 2 && $outdegree >= 2) ||
			($indegree >= 2 && $outdegree == 0)
			);

		# Create intermediate subdivision node
		my $cut = $v . "|cut";
		$g->add_vertex($cut);

		if ($outdegree == 0) {

			# Leaf-side expansion:
			# parents → cut → original leaf
			for my $p (@pred) {
				$g->delete_edge($p, $v);
				$g->add_edge($p, $cut);
			}

			$g->add_edge($cut, $v);

		} else {

			# Internal expansion:
			# original node → cut → children
			for my $c (@succ) {
				$g->delete_edge($v, $c);
				$g->add_edge($cut, $c);
			}

			$g->add_edge($v, $cut);
		}
	}
}
