#######################################################
#################### Run Squirrel ####################
#######################################################

#This scrip runs Squirrel and the boostrap sampling, it does it 10 times (10 different boostrap samples per network)

### This script was written by a human and optimized using Claude LLM

 
base_dir <- path.expand(
  "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final/auxfiles"
)

source_folder <- file.path(base_dir, "sim_out")
python_exe    <- "/Users/006512736/miniconda3/bin/python3"
n_iter        <- 10

# The line number in each python script that sets sim_dir (same in all 4 files)
sim_dir_line <- 453
 
configs <- list(
  list(
    tag                 = "3_rooted_full",
    target_fasta_names  = sprintf("concatalignRealisticSim_400gts_Net%d.fasta", 1:20),
    seqlens_file        = file.path(source_folder, "seq_lens.csv"),
    python_script       = file.path(base_dir, "3SQUIRREL_rooted_full"),
    dir_prefix          = "rootedfull",
    out_py_prefix       = "3SQUIRREL_rooted_full_"
  ),
  list(
    tag                 = "4_rooted_third",
    target_fasta_names  = sprintf("concatalignRealisticSim_400gts_Net%d_third.fasta", 1:20),
    seqlens_file        = file.path(source_folder, "seq_lens_third.csv"),
    python_script       = file.path(base_dir, "4SQUIRREL_rooted_third"),
    dir_prefix          = "rootedthird",
    out_py_prefix       = "4SQUIRREL_rooted_third_"
  ),
  list(
    tag                 = "5_unrooted_full",
    target_fasta_names  = sprintf("concatalignRealisticSim_400gts_Net%d.fasta", 1:20),
    seqlens_file        = file.path(source_folder, "seq_lens.csv"),
    python_script       = file.path(base_dir, "5SQUIRREL_unrooted_full"),
    dir_prefix          = "unrootedfull",
    out_py_prefix       = "5SQUIRREL_unrooted_full_"
  ),
  list(
    tag                 = "6_unrooted_third",
    target_fasta_names  = sprintf("concatalignRealisticSim_400gts_Net%d_third.fasta", 1:20),
    seqlens_file        = file.path(source_folder, "seq_lens_third.csv"),
    python_script       = file.path(base_dir, "6SQUIRREL_unrooted_third"),
    dir_prefix          = "unrootedthird",
    out_py_prefix       = "6SQUIRREL_unrooted_third_"
  )
)

 
for (cfg in configs) {
  fasta_paths <- file.path(source_folder, cfg$target_fasta_names)
  missing_fasta <- fasta_paths[!file.exists(fasta_paths)]
  if (length(missing_fasta) > 0) {
    stop(sprintf(
      "[%s] Missing fasta files:\n%s",
      cfg$tag, paste("  -", missing_fasta, collapse = "\n")
    ))
  }
  if (!file.exists(cfg$seqlens_file)) {
    stop(sprintf("[%s] seqlens file not found: %s", cfg$tag, cfg$seqlens_file))
  }
  if (!file.exists(cfg$python_script)) {
    stop(sprintf("[%s] python script template not found: %s", cfg$tag, cfg$python_script))
  }
}
 
run_one <- function(cfg, i) {
  
  fasta_paths <- file.path(source_folder, cfg$target_fasta_names)
  
  # 1. Create output folder for this config + iteration
  run_dir <- file.path(base_dir, "sim_out", paste0(cfg$dir_prefix, i))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Copy the 20 MSAs and seqlens file into it
  file.copy(fasta_paths, file.path(run_dir, basename(fasta_paths)), overwrite = TRUE)
  file.copy(cfg$seqlens_file, file.path(run_dir, basename(cfg$seqlens_file)), overwrite = TRUE)
  
  # sanity check
  copied_fastas <- list.files(run_dir, pattern = "\\.fasta$")
  if (length(copied_fastas) != length(cfg$target_fasta_names)) {
    stop(sprintf(
      "[%s] iter %d: expected %d fasta files in %s but found %d",
      cfg$tag, i, length(cfg$target_fasta_names), run_dir, length(copied_fastas)
    ))
  }
  
  # 3. Create an edited copy of the python script pointing at this run's folder
  script_lines <- readLines(cfg$python_script, warn = FALSE)
  script_lines[sim_dir_line] <- paste0(
    "sim_dir = base_dir / \"sim_out/", cfg$dir_prefix, i, "\""
  )
  script_copy_path <- file.path(run_dir, paste0(cfg$out_py_prefix, i, ".py"))
  writeLines(script_lines, script_copy_path)
  
  # 4. Run it
  message(sprintf("[%s] Running iteration %d -> %s", cfg$tag, i, run_dir))
  
  result <- system2(
    command = python_exe,
    args    = shQuote(script_copy_path),
    stdout  = file.path(run_dir, "stdout.log"),
    stderr  = file.path(run_dir, "stderr.log")
  )
  
  if (result != 0) {
    warning(sprintf("[%s] iteration %d exited with non-zero status: %d", cfg$tag, i, result))
  }
}

for (i in 6:(n_iter)) {
  for (cfg in configs) {
    run_one(cfg, i)
  }
}

message("All ", n_iter, " interleaved iterations of all 4 pipelines complete.")