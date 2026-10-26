#!/usr/bin/perl 

#HMS-consensus.pl
#2026/04/25
#
# Purpose:
#   Computes the HMS consensus phylogenetic network from Huber, Moulton, and Spillner (2023) from a set of input networks.
#   Input is a CSV-like file with extended Newick strings and optional weights.
#	Output is the uncompressed consensus network as an extended Newick string.
#	To obtain the compressed consensus network as in Huber et al. (2023), contract all edges leaving a reticulation vertex.
#
# Sample usage:
#   perl HMS-consensus.pl --in=input-networks.txt --out=consensus-net.txt


use strict;
use warnings;
use Getopt::Long;	

use Bio::PhyloNetwork;	# For phylogenetic network parsing and manipulation
use Graph;	 # General graph data structure


#######################################################################
# ------------------- COMMAND LINE HANDLING ---------------------------
#######################################################################

# Variables for the command line arguments
my ($input,$output);

# Parse --in and --out options
GetOptions('in=s'=>\$input,'out=s'=>\$output);

# Ensure input file is provided
if (!defined($input)){
    print "Missing network input filename.\n";
    exit(1);
}

# Ensure input file exists
open(my $in, "<", $input) or die "Could not open file '$input'. $!";

# Handle output behavior
our $terminal = 0; 

if (!defined($output)){
	 # No output file → print to terminal
    $terminal = 1;
}

elsif (-e "$output"){
	 # If output exists, ask user what to do
	print "\nThe file $output already exists. \nIf you want to overwrite the file, type OVERWRITE,\nif you want to append the new data to the old file, type APPEND. \nElse the program will die.\n";
	my $answer = <STDIN>;
	chomp($answer);
	
	if($answer eq "OVERWRITE" ){
		unlink ($output) or die;	
	}
	elsif(not $answer eq "APPEND"){
		die;
	}
}



#######################################################################
# ------------------- INPUT PARSING -----------------------------------
#######################################################################

# Read input file of networks
# Format:
#   extended Newick string, optional weight
#   If weight missing → default = 1


open(my $fh, '<', $input) or die "Cannot open $input: $!\n";

my @networks;

# <$fh>;  #discard header line if first line corresponds to "true" network and subsequent lines corresponed to bootstrapped networks

while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g;	# trim whitespace
    next if $line eq '';	# skip empty lines

    my ($newick, $weight);

     # Extract Newick string and optional weight
    if ($line =~ /^(.*?;)\s*(?:,\s*([0-9.eE+-]+)\s*)?$/) {
        $newick = $1;
        $weight = defined $2 ? $2 : 1;
    }
    else {
        die "Line does not look like extended Newick:\n$line\n";
    }

    push @networks, {
        newick => $newick,
        weight => $weight,
    };
	
	#print "Number of networks: ".@networks."\n";
}

close($fh);



#######################################################################
# ------------------- CONSENSUS SET PAIRS COMPUTATION ---------------------------
#######################################################################


my %PAIR;	# Global lookup for S/H pairs (used later)

my %pair_counts;	# Counts weighted occurrences of each S|H pair

# Loop over all networks
foreach my $net (@networks) {
	
    my $weight = $net->{weight} // 1;   # default weight = 1 if not provided
	
	# Extract S/H pairs from network
    my %setpairs = get_SetPairs($net);
	
    my %seen_keys;	  # avoid double counting within one network

    foreach my $node (keys %setpairs) {
		
		 # Sort sets for canonical representation
        my @S_sorted = sort keys %{ $setpairs{$node}{S} };
        my @H_sorted = sort keys %{ $setpairs{$node}{H} };

        my $S_str = join(",", @S_sorted);
        my $H_str = join(",", @H_sorted);
        my $key = "$S_str|$H_str";
		
		next if $seen_keys{$key};

		 # Add weighted contribution
		$pair_counts{$key} += $weight; 
		$seen_keys{$key} = 1;
    }
}

# Compute total weight
my $total_weight = 0;
$total_weight += $_->{weight} // 1 for @networks;



#######################################################################
# ------------------- CONSENSUS FILTERING -----------------------------
#######################################################################

my @consensus_pairs;

foreach my $key (keys %pair_counts) {
    my ($S_str, $H_str) = split(/\|/, $key);
    my $weight_sum = $pair_counts{$key};

    my @S = split(",", $S_str);
    my @H = split(",", $H_str);

    # Majority rules:
    # - If H empty → > 50%
    # - If H non-empty → > 2/3
    if (!@H) {
        if ($weight_sum > $total_weight/2) {
            push @consensus_pairs, { S => \@S, H => \@H, weight => $weight_sum };
        }
    } else {
        if ($weight_sum > (2/3)*$total_weight) {
            push @consensus_pairs, { S => \@S, H => \@H, weight => $weight_sum };
        }
    }
}



#######################################################################
# ------------------- GRAPH (HASSE DIAGRAM) CONSTRUCTION ------------------------------
#######################################################################

# Build full partial order graph
my $g_full = build_graph(\@consensus_pairs, \%PAIR);
   
# Compute Hasse diagram (transitive reduction) 
my $g_hasse  = transitive_reduction($g_full);
   
# Reverse direction (all edged directed away from the root)   
$g_hasse = $g_hasse->transpose_graph;

# Add missing shortcut edges (for details, see Huber et al. 2023)
add_missing_shortcuts($g_hasse,\%PAIR);   

# Clean up graph structure  
relabel_leaf_vertices($g_hasse,\%PAIR); 
uncompress($g_hasse);



#######################################################################
# ------------------- NODE RELABELING ---------------------------------
#######################################################################
   
# Copy graph for renaming
my $g = $g_hasse->deep_copy_graph;

my $i_count = 1; # internal nodes
my $h_count = 1; # reticulation nodes
my %rename;


foreach my $v ($g->vertices) {
	my $out = $g->out_degree($v);
	my $in  = $g->in_degree($v);

	next if $out == 0;	# skip leaves

	# Reticulation node (multiple parents)
	if ($in > 1) {
		$rename{$v} = "#H$h_count";
	    $h_count++;
	} else {
	    $rename{$v} = "I$i_count";
	    $i_count++;
	}
}


# Build renamed graph
my $g2 = Graph->new(directed => 1);

foreach my $v ($g->vertices) {
	my $new = $rename{$v} // $v;  # leaves unchanged
	$g2->add_vertex($new);
}

foreach my $e ($g->edges) {
	my ($src, $dst) = @$e;
	my $new_src = $rename{$src} // $src;
	my $new_dst = $rename{$dst} // $dst;
	$g2->add_edge($new_src, $new_dst);
}



#######################################################################
# ------------------- OUTPUT ------------------------------------------
#######################################################################

# Convert graph back to extended Newick
my $consensusnet = new Bio::PhyloNetwork(-graph => $g2);
my $str=$consensusnet->eNewick();


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
# ------------------- CORE FUNCTION: SET PAIRS ------------------------
#######################################################################  

# Extract S/H sets for each node in a network
sub get_SetPairs {
	
	my ($net) = @_;
	my $newick = $net->{newick};
	die "Newick string undefined!" unless defined $newick;
	
	# Build network object
	our $network = Bio::PhyloNetwork->new(-eNewick => $newick);
	my $graph = $network->graph or die "Graph undefined";

	# Identify leaves
	my @nodes = $graph->vertices;
	my %leaves = map { $_ => 1 } grep { $graph->out_degree($_) == 0 } @nodes;

	# Find root
	my ($root) = grep { $graph->in_degree($_) == 0 } @nodes;
	die "Root not found!" unless defined $root;

	 # Compute clusters (descendant leaves)
    my %Cluster;
    foreach my $u (@nodes) {
        my %cluster;
        $cluster{$u} = 1 if exists $leaves{$u};
		
        foreach my $v ($graph->all_successors($u)) {
            $cluster{$v} = 1 if exists $leaves{$v};
        }
		
        $Cluster{$u} = \%cluster;
	}



	# Helper: enumerate all root→leaf paths
	sub all_paths {
		my ($graph, $current, $target, $path) = @_;
		$path ||= [];
	
		my @paths;
		my @new_path = (@$path, $current);
	
		if ($current eq $target) {
			push @paths, [@new_path];
			return @paths;
		}

		foreach my $child ($graph->successors($current)) {
			push @paths, all_paths($graph, $child, $target, \@new_path);
		}
		return @paths;
	}

	# Compute S and H sets:
	    # S = leaves where all paths go through node
	    # H = leaves where some paths bypass node
	my %setpairs;
	
	foreach my $u (@nodes) {
		
		my %S;
		my %H;

		foreach my $leaf (keys %{ $Cluster{$u} }) {
			
			my @paths = all_paths($graph, $root, $leaf);
			
			my $all_through_u = 1;
			
			foreach my $path (@paths) {
				unless (grep { $_ eq $u } @$path) {
					$all_through_u = 0;
					last;
				}
			}

			if ($all_through_u) {
				$S{$leaf} = 1;
			} else {
				$H{$leaf} = 1;
			}
		}

		$setpairs{$u} = { S => \%S, H => \%H };
	}

	return %setpairs;
}



#######################################################################
# ------------------- BASIC SET OPERATIONS ----------------------------
#######################################################################

# Return union of two sets (represented as hashrefs)
sub set_union {
	my ($A,$B) = @_;
	my %u = (%$A, %$B);  # merge keys
	return \%u;
}


# Check if A ⊆ B
sub is_subset {
	my ($A,$B) = @_;
	foreach my $x (keys %$A) {
		return 0 unless exists $B->{$x};
	}
	return 1;
}


# Check if A ⊂ B (proper subset)
sub is_proper_subset {
	my ($A,$B) = @_;
	return 0 unless is_subset($A,$B);
	return scalar(keys %$A) < scalar(keys %$B);
}


# Check if two sets are equal
sub set_equal {
	my ($A, $B) = @_;

	return 0 unless keys(%$A) == keys(%$B);   # same size

	for my $k (keys %$A) {
		return 0 unless exists $B->{$k};      # every element matches
	}

	return 1;
}



#######################################################################
# ------------------- ORDER RELATION ON SET PAIRS -------------------------
#######################################################################

# Define partial order between two (S,H) pairs
# Returns 1 if p1 < p2
#
# This is the key relation used to build the Hasse diagram.
sub pair_less {
	my ($p1,$p2) = @_;#

	# Convert arrays → set representation
	my %S1 = map { $_ => 1 } @{$p1->{S}};
	my %H1 = map { $_ => 1 } @{$p1->{H}};
	my %S2 = map { $_ => 1 } @{$p2->{S}};
	my %H2 = map { $_ => 1 } @{$p2->{H}};#

	# U1 = S1 ∪ H1
	my $U1 = set_union(\%S1,\%H1);###

	# Condition (a): U1 ⊆ S2
	return 1 if is_subset($U1,\%S2);
	
	# Condition (b): U1 ⊆ H2
	return 1 if is_subset($U1,\%H2);
	
	# Condition (c):
	# Same H-set AND S1 ⊂ S2
	if (set_equal(\%H1,\%H2) && is_proper_subset(\%S1,\%S2)) {
		return 1;
    }

	return 0;
}



#######################################################################
# ------------------- SERIALIZATION -----------------------------------
#######################################################################

# Convert a pair {S,H} into a canonical string "S|H"
# Used as graph vertex identifier
sub serialize_pair {
	my ($p) = @_;
	return $p if !ref($p);	 # already a string

	my $S = join(",", sort @{$p->{S}});
	my $H = join(",", sort @{$p->{H}});

	$S //= '';
	$H //= '';

	return "$S|$H";
}



#######################################################################
# ------------------- GRAPH (HASSE DIAGRAM) CONSTRUCTION ------------------------------
#######################################################################

# Build directed graph encoding partial order between pairs
sub build_graph {
	my ($pairs_ref, $PAIR) = @_;   # receive references correctly

	my $g = Graph->new(directed => 1);

	# ---- Add vertices ----
    foreach my $p (@$pairs_ref) {

        my $key = serialize_pair($p);
        $g->add_vertex($key);

		# Store original S/H sets for later lookup
		$PAIR->{$key} = {
			S => $p->{S},
			H => $p->{H},
		};
	}

	  # ---- Add edges according to pair_less relation ----
	  for my $i (0 .. $#$pairs_ref) {
		for my $j (0 .. $#$pairs_ref) {
			next if $i == $j;

			my $pi = $pairs_ref->[$i];
			my $pj = $pairs_ref->[$j];

			if (pair_less($pi, $pj)) {
				my $u = serialize_pair($pi);
				my $v = serialize_pair($pj);
				$g->add_edge($u, $v);
			}
		}
	}

	return $g;
}



#######################################################################
# ------------------- TRANSITIVE REDUCTION ----------------------------
#######################################################################

# Check if target reachable from start WITHOUT using edge (skip_u → skip_v)
sub reachable_without_edge {
	my ($g,$start,$target,$skip_u,$skip_v) = @_;

	my %seen;
	my @stack = ($start);

	while (@stack) {
		my $x = pop @stack;
		next if $seen{$x}++;

		foreach my $y ($g->successors($x)) {

			# Ignore the tested edge
			next if ($x eq $skip_u && $y eq $skip_v);

			return 1 if $y eq $target;
			push @stack, $y;
		}
	}
	return 0;
}



# Compute Hasse diagram (remove transitive edges)
sub transitive_reduction {
	my ($g) = @_;

	my $reduced = $g->copy;

	foreach my $e ($g->edges) {
		my ($u,$v) = @$e;

		# If alternate path exists, edge is redundant
		if (reachable_without_edge($g,$u,$v,$u,$v)) {
			$reduced->delete_edge($u,$v);
		}
	}

	return $reduced;
}



#######################################################################
# ------------------- SHORTCUT EDGES ----------------------------------
#######################################################################

# Helper
# Canonical string representation of a set
sub set_key {
	my ($aref) = @_;
	return join(",", sort @$aref); 
}


# Add edges required by the consensus construction that are
# not present in the pure Hasse diagram
sub add_missing_shortcuts {
	
	my ($g,$PAIR) = @_;

	# Work on frozen copy to avoid modifying during traversal
	my $base = $g->copy;

	foreach my $v ($base->vertices) {
		
		# Get parent(s)
		my @parents = $base->predecessors($v);
		next unless @parents == 1;

		my $p = $parents[0];

		# Condition: H(parent) = S(child)
		my $S_v = set_key($PAIR->{$v}{S});
		my $H_p = set_key($PAIR->{$p}{H});

		next unless $H_p eq $S_v;
		
		# Walk upward as long as H stays constant
		my $curr = $p;

		while (1) {
			my @pp = $base->predecessors($curr);
			last unless @pp == 1;

			my $next = $pp[0];

			my $H_curr = set_key($PAIR->{$curr}{H});
			my $H_next = set_key($PAIR->{$next}{H});

			last unless $H_curr eq $H_p;

			$curr = $next;
		}

		 # Add shortcut edge
		 $g->add_edge($curr,$v) unless $g->has_edge($curr,$v);
    }
}



#######################################################################
# ------------------- GRAPH CLEANUP -----------------------------------
#######################################################################


# Insert subdivision nodes ("cut" vertices) to return an uncompressed network
sub uncompress {
	my ($g) = @_;

	my @vertices = $g->vertices;

	for my $v (@vertices) {

		my @pred = $g->predecessors($v);
		my @succ = $g->successors($v);

		my $indegree  = scalar @pred;
		my $outdegree = scalar @succ;

		# Only process specific structural cases
		next unless (
			($indegree >= 2 && $outdegree >= 2) ||
			($indegree >= 2 && $outdegree == 0)
		);

		# Create subdivision node
		my $cut = $v . "|cut";
		$g->add_vertex($cut);

		if ($outdegree == 0){
			# Leaf case: parents → cut → v
			for my $p (@pred) {
				$g->delete_edge($p,$v);
				$g->add_edge($p,$cut);
			}
			$g->add_edge($cut,$v);

		} else {
			# Internal case: v → cut → children
			for my $c (@succ) {
				$g->delete_edge($v,$c);
				$g->add_edge($cut,$c);
			}
			$g->add_edge($v,$cut);
		}
	}
}



#######################################################################
# ------------------- LEAF RELABELING ---------------------------------
#######################################################################

# Replace encoded leaf nodes "x|" with plain taxon labels
sub relabel_leaf_vertices {
	my ($g,$PAIR) = @_;
	
	my @vertices = $g->vertices;
	
	foreach my $v (@vertices) {
		
		my ($S,$H) = split(/\|/,$v,-1);
		
		my @S_elems = grep { $_ ne '' } split(/,/,$S);
		my @H_elems = grep { $_ ne '' } split(/,/,$H);

		# Detect leaf: S = {x}, H = ∅
		next unless @S_elems == 1 && @H_elems == 0;
		
		my $leaf_name = $S_elems[0];

		next if $v eq $leaf_name;

		# Replace vertex name
		$g->add_vertex($leaf_name);

		# Redirect edges
		foreach my $p ($g->predecessors($v)) {
			$g->add_edge($p,$leaf_name);
		}
		foreach my $c ($g->successors($v)) {
			$g->add_edge($leaf_name,$c);
		}

		$g->delete_vertex($v);

		# Update lookup table
		$PAIR->{$leaf_name} = $PAIR->{$v};
		delete $PAIR->{$v};
	}
}
