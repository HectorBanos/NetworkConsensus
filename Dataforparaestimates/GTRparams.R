library(castor)
library(ape)
library(ggplot2)
library(phangorn)
library(adephylo)


# =============================================================================
# Working directory
# =============================================================================

setwd("~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/iqtree-3.0.1-macOS/bin")
 
okabe_ito <- c(
  "#E69F00", 
  "#56B4E9",  
  "#009E73", 
  "#F0E442",  
  "#0072B2",  
  "#D55E00",  
  "#CC79A7",  
  "#000000", 
  "#999999",  
  "#88CCEE",  
  "#44AA99",  
  "#117733",  
  "#332288",  
  "#DDCC77",  
  "#CC6677", 
  "#AA4499"  
)
 
getalpha <- function(file) {
  
  filename <- paste0(file, ".iqtree")
  lines <- readLines(filename)
  
  start <- which(
    grepl("Model of rate heterogeneity:", lines, fixed = TRUE)
  ) + 1
  
  alpha <- as.numeric(strsplit(lines[start], "\\s+")[[1]][4])
  
  return(alpha)
}

# Extract six GTR exchangeability parameters
getrates <- function(file) {
  
  filename <- paste0(file, ".iqtree")
  lines <- readLines(filename)
  
  start <- which(
    grepl("Rate parameter R:", lines, fixed = TRUE)
  ) + 2
  
  end <- which(
    grepl(
      "State frequencies: (empirical counts from alignment)",
      lines,
      fixed = TRUE
    )
  ) - 2
  
  rates <- numeric()
  
  counter <- 1
  for (i in start:end) {
    
    rates[counter] <- as.numeric(
      gsub("=", "", substr(lines[i],
                           nchar(lines[i]) - 5,
                           nchar(lines[i])))
    )
    
    counter <- counter + 1
  }
  
  return(rates)
}

# Extract nucleotide frequencies (A,C,G,T)
getpi <- function(file) {
  
  filename <- paste0(file, ".iqtree")
  lines <- readLines(filename)
  
  start <- which(
    grepl(
      "State frequencies: (empirical counts from alignment)",
      lines,
      fixed = TRUE
    )
  ) + 2
  
  end <- which(
    grepl("Rate matrix Q:", lines, fixed = TRUE)
  ) - 2
  
  pie <- numeric()
  
  counter <- 1
  for (i in start:end) {
    
    pie[counter] <- as.numeric(
      gsub("=", "", substr(lines[i],
                           nchar(lines[i]) - 5,
                           nchar(lines[i])))
    )
    
    counter <- counter + 1
  }
  
  return(pie)
}


 
genes <- c(
  "10PrimRAG1",
  "10PrimRAG2",
  "10PrimIRBP",
  "10PrimBRCA1",
  "10PrimAPOB",
  "10PrimVWF",
  "10PrimGHR",
  "10PrimMC1R",
  "10PrimTYR",
  "10PrimMITF"
)

 
atrees <- numeric(length(genes))

pitrees <- matrix(
  0,
  nrow = length(genes),
  ncol = 4
)

ratetrees <- matrix(
  0,
  nrow = length(genes),
  ncol = 6
)

for (i in seq_along(genes)) {
  
  file <- paste0(genes[i], "v2.fa")
  
  atrees[i]      <- getalpha(file)
  pitrees[i, ]   <- getpi(file)
  ratetrees[i, ] <- getrates(file)
}


gene_labels <- gsub("10Prim", "", genes)

gene_shapes <- c(
  16, 17, 15, 18, 8,
  3, 4, 7, 9, 10
)

common_theme <- theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
# =============================================================================
# Figure 1: GTR exchangeability parameters
#
# Parameters correspond to:
# AC, AG, AT, CG, CT, GT
# =============================================================================

exchange_df <- data.frame(
  genes = factor(
    rep(seq_len(nrow(ratetrees)), times = ncol(ratetrees)),
    labels = gene_labels
  ),
  parameter = factor(
    rep(seq_len(ncol(ratetrees)), each = nrow(ratetrees))
  ),
  value = as.vector(ratetrees)
)

ggplot(exchange_df, aes(x = parameter, y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey85",
    color = "black"
  ) +
  geom_point(
    aes(color = genes, shape = genes),
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.08)
  ) +
  scale_colour_manual(
    values = okabe_ito,
    name = "Gene"
  ) +
  scale_shape_manual(
    values = gene_shapes,
    name = "Gene"
  ) +
  labs(
    title = "Exchangeability across primate genes",
    x = "",
    y = "Exchangeability rate"
  ) +
  common_theme

median_exch <- aggregate(value ~ parameter, exchange_df, median)
mean_exch   <- aggregate(value ~ parameter, exchange_df, mean)
median_exch 

# =============================================================================
# Figure 2: Nucleotide frequencies
#
# Frequencies correspond to:
# pi_A, pi_C, pi_G, pi_T
# =============================================================================

basefreq_df <- data.frame(
  genes = factor(
    rep(seq_len(nrow(pitrees)), times = ncol(pitrees)),
    labels = gene_labels
  ),
  nucleotide = factor(
    rep(seq_len(ncol(pitrees)), each = nrow(pitrees))
  ),
  value = as.vector(pitrees)
)

nucleotide_labels <- c(
  expression(pi[A]),
  expression(pi[C]),
  expression(pi[G]),
  expression(pi[T])
)

ggplot(basefreq_df, aes(x = nucleotide, y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey85",
    color = "black"
  ) +
  geom_point(
    aes(color = genes, shape = genes),
    size = 1.8,
    alpha = 0.7
  ) +
  scale_colour_manual(
    values = okabe_ito,
    name = "Gene"
  ) +
  scale_shape_manual(
    values = gene_shapes,
    name = "Gene"
  ) +
  scale_x_discrete(labels = nucleotide_labels) +
  labs(
    title = "Base frequencies across primate genes",
    x = "Nucleotide",
    y = "Frequency"
  ) +
  common_theme

normalized_median_basedist <- with(
  aggregate(value ~ nucleotide, basefreq_df, median),
  value / sum(value)
)

normalized_mean_basedist <- with(
  aggregate(value ~ nucleotide, basefreq_df, mean),
  value / sum(value)
)

normalized_median_basedist
# =============================================================================
# Figure 3: Gamma shape parameter (alpha)
#
# Smaller alpha values indicate stronger among-site rate heterogeneity.
# =============================================================================

alpha_df <- data.frame(
  genes = factor(gene_labels),
  group = factor(rep("All genes", length(atrees))),
  value = atrees
)

ggplot(alpha_df, aes(x = group, y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey85",
    color = "black"
  ) +
  geom_point(
    aes(color = genes, shape = genes),
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.08)
  ) +
  scale_colour_manual(
    values = okabe_ito,
    name = "Gene"
  ) +
  scale_shape_manual(
    values = gene_shapes,
    name = "Gene"
  ) +
  labs(
    title = "Gamma shape parameter across primate genes",
    x = "",
    y = expression(alpha)
  ) +
  common_theme

median_rates <- aggregate(value ~ group, alpha_df, median)
mean_rates   <- aggregate(value ~ group, alpha_df, mean)
 
median_rates

############## Papionini data used for getting the Diameters and the MaxMin Ratio


setwd("~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline/Dataforparaestimates")

# =============================================================================
# Target taxa (ingroup definition)
# =============================================================================

taxa <- c(
  "Macaca_nemestrina",
  "Macaca_fascicularis",
  "Macaca_mulatta",
  "Mandrillus_leucophaeus",
  "Cercocebus_atys",
  "Papio_anubis",
  "Theropithecus_gelada",
  "Callithrix_jacchus"
)

# =============================================================================
# Load trees
# =============================================================================

papio <- read.tree("papionini.treefile")

# Storage vectors
max_ratio_papio <- numeric(length(papio))
tree_diameters  <- numeric(length(papio))

# =============================================================================
# Loop over trees
# =============================================================================

for (i in seq_along(papio)) {
  
  # ---------------------------------------------------------------------------
  # Step 1: keep only ingroup taxa
  # ---------------------------------------------------------------------------
  
  pruned_tree <- drop.tip(
    papio[[i]],
    tip = papio[[i]]$tip.label[!papio[[i]]$tip.label %in% taxa]
  )
  
  # ---------------------------------------------------------------------------
  # Step 2: midpoint root the tree
  # ---------------------------------------------------------------------------
  
  ingroup_tree <- midpoint(pruned_tree)
  
  # ---------------------------------------------------------------------------
  # Step 3: compute tree diameter (max pairwise distance)
  # ---------------------------------------------------------------------------
  
  tree_diameters[i] <- max(cophenetic(ingroup_tree))
  
  # ---------------------------------------------------------------------------
  # Step 4: root-to-tip distance ratio (imbalance proxy)
  # ---------------------------------------------------------------------------
  
  root_distances <- get_all_distances_to_root(ingroup_tree)[
    seq_along(ingroup_tree$tip.label)
  ]
  
  max_ratio_papio[i] <- max(root_distances) / min(root_distances)
}

# =============================================================================
# Filter extreme tree diameters
# =============================================================================

sum(tree_diameters > 0.14)

tree_diameters <- tree_diameters[tree_diameters < 0.14]

# Note: ~140 trees removed (out of 1730) as extreme values

write.csv(
  x = tree_diameters,
  file = "papionini_diameters.csv",
  row.names = FALSE
)

# =============================================================================
# Summary statistic
# =============================================================================

median(max_ratio_papio)
 
 
