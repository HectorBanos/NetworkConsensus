
########################################################################### 
#################### Create excel files values  ###########################
########################################################################### 

### This script was written by a human and optimized using Claude LLM



library(openxlsx)

base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-final"
sim_out  <- file.path(base_dir, "sim_out")
n_reps   <- 10

 
blob_header1_simple <- c(
  "True", "", "Blob", "", "Circular ord", "", "Hybrid", "", "", "SQUIRREL",
  "", "Blob", "","", "Circular ord", "","", "Hybrid", "",  "", "",   "Extended","",
  "Blob", "", "", "Circular ord", "", "", "Hybrid", "",""
)
blob_header2_simple <- c(
  "Network", "Blob size", "TS", "BP", "TS", "BP", "TS", "BP",
  "",
  "Network", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Net",
  "TS", "BP", "TS vs True Tob",
  "",
  "Network 1", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob"
)

 
blob_header1_extended <- c(
  "True", "", "Blob", "", "Circular ord", "", "Hybrid", "", "", "SQUIRREL",
  "",
  "Blob", "", "", "Circular ord", "", "", "Hybrid", "", "", "", "Extended",
  "",
  "Blob", "", "", "Circular ord", "", "", "Hybrid", "",
  "", "",
  "Huber",
  "", "",
  "Blob", "", "", "Circular ord", "", "", "Hybrid", "", "",  "HWC",
  "",
  "Blob", "", "", "Circular ord", "", "", "Hybrid", ""
)
blob_header2_extended <- c(
  "Network", "Blob size", "TS", "BP", "TS", "BP", "TS", "BP",
  "",
  "Network", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Net",
  "TS", "BP", "TS vs True Tob",
  "",
  "Network", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "",
  "Network", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "",
  "Network", "Blob size", "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob",
  "TS", "BP", "TS vs True Tob"
)

 
edge_header1_extended <- c(
  "",
  "# of branches", "", "", "", "", "",
  "RF dist vs True Tree", "", "", "", "",
  "Extended", "", "",
  "Squirrel", "", "",
  "Huber", "", "",
  "Hardwire", "", "",
  "Majority", "", "",
  "AVG RF to bootstrap", "", "", "", "", ""
)
edge_header2_extended <- c(
  "Network",
  "TRUE", "Extended", "Squirrel", "Huber",  "Hardwire","Majority",
  "Extended", "Squirrel", "Huber", "Majority", "Hardwire",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "TRUE", "Extended", "Squirrel", "Huber",  "Hardwire","Majority"
)

 
edge_header1_simple <- c(
  "",
  "# of branches", "", "", "",
  "RF dist vs True Tree", "", "",
  "Extended", "","",
  "Squirrel", "","",
  "Majority", "","",
  "AVG RF to bootstrap", "", "", ""
)
edge_header2_simple <- c(
  "Network",
  "TRUE", "Extended", "Squirrel", "Majority",
  "Extended", "Squirrel", "Majority",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "TRUE", "Extended", "Squirrel", "Majority"
)

parse_canon <- function(canon_str) {
  blocks <- strsplit(canon_str, " / ", fixed = TRUE)[[1]]
  lapply(blocks, function(b) strsplit(b, "*", fixed = TRUE)[[1]])
}


restrict_partition <- function(blocks, taxa_keep) {
  blocks <- lapply(blocks, function(b) sort(intersect(b, taxa_keep)))
  blocks <- blocks[lengths(blocks) > 0]
  blocks[order(vapply(blocks, `[`, character(1), 1))]
}


canon_compatible <- function(canon1, canon2) {
  b1 <- parse_canon(canon1)
  b2 <- parse_canon(canon2)
  common <- intersect(unlist(b1), unlist(b2))
  if (length(common) == 0) return(FALSE)
  identical(restrict_partition(b1, common), restrict_partition(b2, common))
}


merge_zero_rows <- function(df) {
  score_cols <- setdiff(names(df), c("canon", "n_partitions"))
  first_col <- score_cols[1]
  
  repeat {
    zero_rows <- which(df[[first_col]] == 0)
    if (length(zero_rows) == 0) break
    
    i <- zero_rows[which.min(df$n_partitions[zero_rows])]
    s <- df$n_partitions[i]
    
    nonzero_cols_i <- score_cols[as.numeric(df[i, score_cols]) != 0]
    if (length(nonzero_cols_i) == 0) break
    
    candidates <- which(
      df$n_partitions > s &
        apply(df[, nonzero_cols_i, drop = FALSE], 1,
              function(r) all(as.numeric(r) == 0)) &
        vapply(df$canon, canon_compatible, logical(1), canon2 = df$canon[i])
    )
    
    if (length(candidates) == 0) break
    
    j <- candidates[which.min(df$n_partitions[candidates])]
    
    for (col in nonzero_cols_i) {
      df[j, col] <- df[i, col]
    }
    df <- df[-i, ]
  }
  
  rownames(df) <- NULL
  return(df)
}


write_blob_workbook <- function(L, header1, header2, outfile) {
  
  wb <- createWorkbook()
  addWorksheet(wb, "Matrices")
  
  start_row   <- 4
  header_rows <- 2
  
  for (i in seq_along(L)) {
    
    writeData(wb, "Matrices", t(as.matrix(header1)),
              startRow = start_row, startCol = 1, colNames = FALSE)
    
    writeData(wb, "Matrices", t(as.matrix(header2)),
              startRow = start_row + 1, startCol = 1, colNames = FALSE)
    
    matrix_row <- start_row + header_rows
    start_col  <- 1
    max_rows   <- 0
    
    for (j in seq_along(L[[i]])) {
      
     
      L[[i]][[j]][is.na(L[[i]][[j]])] <- 0
      M <- L[[i]][[j]]
       
      if (is.data.frame(M) && all(c("canon", "n_partitions") %in% names(M))) {
        score_cols <- setdiff(names(M), c("canon", "n_partitions"))
        if (length(score_cols) > 0 && any(M[[score_cols[1]]] == 0, na.rm = TRUE)) {
          M <- merge_zero_rows(M)
        }
      }
      
      M[,1] <- 1:length(M[,1])
      M <- apply(M, 2, as.numeric)
      M <- as.matrix(M)
      if (dim(M)[2] == 1) { M <- t(M) }
      
      writeData(
        wb, "Matrices", as.data.frame(M),
        startRow = matrix_row, startCol = start_col,
        colNames = FALSE, rowNames = FALSE
      )
      
      max_rows  <- max(max_rows, nrow(M))
      start_col <- start_col + ncol(M) + 1
    }
    
    start_row <- matrix_row + max_rows + 4
  }
  
  saveWorkbook(wb, outfile, overwrite = TRUE)
}

write_edge_workbook <- function(M, header1, header2, outfile) {
  
  M <- cbind((1:(dim(M)[1])), M)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  
  writeData(wb, "Results", t(as.matrix(header1)),
            startRow = 1, startCol = 1, colNames = FALSE)
  
  writeData(wb, "Results", t(as.matrix(header2)),
            startRow = 2, startCol = 1, colNames = FALSE)
  
  writeData(
    wb, "Results", as.data.frame(M),
    startRow = 3, startCol = 1,
    colNames = FALSE, rowNames = FALSE
  )
  
  header_style <- createStyle(textDecoration = "bold", halign = "center")
  addStyle(
    wb, "Results", header_style,
    rows = 1:2, cols = 1:max(length(header1), ncol(M)),
    gridExpand = TRUE, stack = TRUE
  )
  
  saveWorkbook(wb, outfile, overwrite = TRUE)
}

 
for (rep_i in 1:n_reps) {
  
  message("=== Excel export: replicate ", rep_i, " ===")
  
   
  dir_rf <- file.path(sim_out, paste0("rootedfull", rep_i))
  load(file.path(dir_rf, "BlobAnalysis_full_outgroup.Rdata")); result_rf <- result
  load(file.path(dir_rf, "EdgeAnalysis_full_outgroup.Rdata"))   
  
  write_blob_workbook(result_rf, blob_header1_extended, blob_header2_extended,
                      file.path(dir_rf, "BlobAnalysis_full_outgroup.xlsx"))
  write_edge_workbook(MM1, edge_header1_extended, edge_header2_extended,
                      file.path(dir_rf, "EdgeAnalysis_full_outgroup.xlsx"))
  
   
  dir_rt <- file.path(sim_out, paste0("rootedthird", rep_i))
  load(file.path(dir_rt, "BlobAnalysis_third_outgroup.Rdata")); result_rt <- result
  load(file.path(dir_rt, "EdgeAnalysis_third_outgroup.Rdata"))   
  
  write_blob_workbook(result_rt, blob_header1_extended, blob_header2_extended,
                      file.path(dir_rt, "BlobAnalysis_third_outgroup.xlsx"))
  write_edge_workbook(MM2, edge_header1_extended, edge_header2_extended,
                      file.path(dir_rt, "EdgeAnalysis_third_outgroup.xlsx"))
  
   
  dir_uf <- file.path(sim_out, paste0("unrootedfull", rep_i))
  load(file.path(dir_uf, "BlobAnalysis_full_nooutgroup.Rdata")); result_uf <- result
  load(file.path(dir_uf, "EdgeAnalysis_full_nooutgroup.Rdata"))   
  
  write_blob_workbook(result_uf, blob_header1_simple, blob_header2_simple,
                      file.path(dir_uf, "BlobAnalysis_full_nooutgroup.xlsx"))
  write_edge_workbook(MM3, edge_header1_simple, edge_header2_simple,
                      file.path(dir_uf, "EdgeAnalysis_full_nooutgroup.xlsx"))
  
   
  dir_ut <- file.path(sim_out, paste0("unrootedthird", rep_i))
  load(file.path(dir_ut, "BlobAnalysis_third_nooutgroup.Rdata")); result_ut <- result
  load(file.path(dir_ut, "EdgeAnalysis_third_nooutgroup.Rdata"))  # loads MM4
  
  write_blob_workbook(result_ut, blob_header1_simple, blob_header2_simple,
                      file.path(dir_ut, "BlobAnalysis_third_nooutgroup.xlsx"))
  write_edge_workbook(MM4, edge_header1_simple, edge_header2_simple,
                      file.path(dir_ut, "EdgeAnalysis_third_nooutgroup.xlsx"))
}

message("Done. Wrote 2 xlsx files into each of the 4 folders for all ", n_reps, " replicates.")