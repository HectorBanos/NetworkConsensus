library(openxlsx)
base_dir <- "~/Library/CloudStorage/OneDrive-CaliforniaStateUniversitySanBernardino/Documents/ToBolivier/R/Pipeline/CleanPipeline-thirdv2"
location <- file.path(base_dir, "sim_out")
setwd(location)
file1=paste0(location,'/BlobAnalysis_full_nooutgroup.Rdata')
file2=paste0(location,'/BlobAnalysis_third_nooutgroup.Rdata')

file3=paste0(location,'/BlobAnalysis_full_outgroup.Rdata')
file4=paste0(location,'/BlobAnalysis_third_outgroup.Rdata')

filenames=c('BlobAnalysis_full_nooutgroup','BlobAnalysis_third_nooutgroup','BlobAnalysis_full_outgroup','BlobAnalysis_third_outgroup')

load(file1)
result1=result
load(file2)
result2=result
load(file3)
result3=result
load(file4)
result4=result


header1 <- c(
  "", "", "Blob", "", "Circular ord", "", "Hybrid", "", "", "SQUIRREL",
  "", "Blob", "","", "Circular ord", "","", "Hybrid", "",  "", "",   "Extended","", 
  "Blob", "", "", "Circular ord", "", "", "Hybrid", "",""
)

header2 <- c(
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
 
for(kk in 1:2)
{
    L=list(result1,result2,result3,result4)[[kk]]
    wb <- createWorkbook()
    addWorksheet(wb, "Matrices")
    start_row <- 4
    header_rows <- 2
    header_rows <- 2
    
    for (i in seq_along(L)) {
      
      ## write headers
      writeData(wb, "Matrices",
                t(as.matrix(header1)),
                startRow = start_row,
                startCol = 1,
                colNames = FALSE)
      
      writeData(wb, "Matrices",
                t(as.matrix(header2)),
                startRow = start_row + 1,
                startCol = 1,
                colNames = FALSE)
      
      ## matrices start below the headers
      matrix_row <- start_row + header_rows
      
      start_col <- 1
      max_rows <- 0
      
      for (j in seq_along(L[[i]])) {
        
        M <- L[[i]][[j]]
        M[,1] <- 1:length(M[,1]) 
        M <- apply(M, 2, as.numeric)
 
        M <- as.matrix(M)
        if(dim(M)[2]==1){M <- t(M)}
        
        writeData(
          wb,
          "Matrices",
          as.data.frame(M),
          startRow = matrix_row,
          startCol = start_col,
          colNames = FALSE,
          rowNames = FALSE
        )
        
        max_rows <- max(max_rows, nrow(M))
        
        start_col <- start_col + ncol(M) + 1
      }
      
      start_row <- matrix_row + max_rows + 4
    }
      saveWorkbook(wb, paste0(filenames[kk],".xlsx"), overwrite = TRUE)
}
















header1 <- c(
  "", "", "Blob", "", "Circular ord", "", "Hybrid", "", "", "SQUIRREL",
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

header2 <- c(
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
 
for (kk in 3:4) {
  
  L <- list(result1, result2, result3, result4)[[kk]]
  
  wb <- createWorkbook()
  addWorksheet(wb, "Matrices")
  
  start_row <- 4
  header_rows <- 2
  
  for (i in seq_along(L)) {
 
    writeData(
      wb,
      "Matrices",
      t(as.matrix(header1)),
      startRow = start_row,
      startCol = 1,
      colNames = FALSE
    )
    
    writeData(
      wb,
      "Matrices",
      t(as.matrix(header2)),
      startRow = start_row + 1,
      startCol = 1,
      colNames = FALSE
    )
  
    matrix_row <- start_row + header_rows
    
    start_col <- 1
    max_rows <- 0
 
    for (j in seq_along(L[[i]])) {
 
      M <- L[[i]][[j]]
      M[,1] <- 1:length(M[,1]) 
      M <- apply(M, 2, as.numeric)
       M <- as.matrix(M)
       if(dim(M)[2]==1){M <- t(M)}
       
      writeData(
        wb,
        "Matrices",
        as.data.frame(M),
        startRow = matrix_row,
        startCol = start_col,
        colNames = FALSE,
        rowNames = FALSE
      )
      
      max_rows <- max(max_rows, nrow(M))
      start_col <- start_col + ncol(M) + 1
    }
 
    start_row <- matrix_row + max_rows + 4
  }
 
  saveWorkbook(
    wb,
    paste0(filenames[kk], ".xlsx"),
    overwrite = TRUE
  )
}










file1=paste0(location,'/EdgeAnalysis_full_nooutgroup.Rdata')
file2=paste0(location,'/EdgeAnalysis_third_nooutgroup.Rdata')

file3=paste0(location,'/EdgeAnalysis_full_outgroup.Rdata')
file4=paste0(location,'/EdgeAnalysis_third_outgroup.Rdata')

filenames=c('EdgeAnalysis_full_nooutgroup','EdgeAnalysis_third_nooutgroup','EdgeAnalysis_full_outgroup','EdgeAnalysis_third_outgroup')

load(file1)
 
load(file2)
 
load(file3)
 
load(file4)
 



header1 <- c(
    "",
  "# of branches", "", "", "", "", "",
  "RF dist vs True Tree", "", "", "", "",
  "Extended", "", "",
  "Squirrel", "", "",
  "Huber", "", "",
  "Majority", "", "",
  "Hardwire", "", "",
  "AVG RF to bootstrap", "", "", "", "", ""
)

header2 <- c(
  
  "Network",
  "TRUE", "Extended", "Squirrel", "Huber", "Majority", "Hardwire",
  "Extended", "Squirrel", "Huber", "Majority", "Hardwire",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "TRUE", "Extended", "Squirrel", "Huber", "Majority", "Hardwire"
)

for( kk in 3:4)
{
  M <- list(MM3, MM4, MM1, MM2)[[kk]] 
  M=cbind((1:(dim(M)[1])),M)
  
 
    wb <- createWorkbook()
    addWorksheet(wb, "Results")
     
    writeData(
      wb, "Results",
      t(as.matrix(header1)),
      startRow = 1,
      startCol = 1,
      colNames = FALSE
    )
    
    writeData(
      wb, "Results",
      t(as.matrix(header2)),
      startRow = 2,
      startCol = 1,
      colNames = FALSE
    )
    
     
    writeData(
      wb,
      "Results",
      as.data.frame(M),
      startRow = 3,
      startCol = 1,
      colNames = FALSE,
      rowNames = FALSE
    )
    
     
    header_style <- createStyle(
      textDecoration = "bold",
      halign = "center"
    )
    
    addStyle(
      wb,
      "Results",
      header_style,
      rows = 1:2,
      cols = 1:max(length(header1), ncol(M)),
      gridExpand = TRUE,
      stack = TRUE
    )
     
    saveWorkbook(
      wb,
      paste0(filenames[kk], ".xlsx"),
      overwrite = TRUE
    )
}









header1 <- c(
   "",
  "# of branches", "", "", "",
  "RF dist vs True Tree", "", "",
  "Extended", "","",
  "Squirrel", "","",
  "Majority", "","",
  "AVG RF to bootstrap", "", "", ""
)
header2 <- c(
  
  "Network",
  "TRUE", "Extended", "Squirrel", "Majority",
  "Extended", "Squirrel", "Majority",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "Matching Splits", "T2 error", "T1 error",
  "TRUE", "Extended", "Squirrel", "Majority"
)

for( kk in 1:2)
{
  M <- list(MM3, MM4, MM1, MM2)[[kk]] 
  M=cbind((1:(dim(M)[1])),M)
   
 
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  
  writeData(
    wb, "Results",
    t(as.matrix(header1)),
    startRow = 1,
    startCol = 1,
    colNames = FALSE
  )
  
  writeData(
    wb, "Results",
    t(as.matrix(header2)),
    startRow = 2,
    startCol = 1,
    colNames = FALSE
  )
  
  
  writeData(
    wb,
    "Results",
    as.data.frame(M),
    startRow = 3,
    startCol = 1,
    colNames = FALSE,
    rowNames = FALSE
  )
  
  
  header_style <- createStyle(
    textDecoration = "bold",
    halign = "center"
  )
  
  addStyle(
    wb,
    "Results",
    header_style,
    rows = 1:2,
    cols = 1:max(length(header1), ncol(M)),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  saveWorkbook(
    wb,
    paste0(filenames[kk], ".xlsx"),
    overwrite = TRUE
  )
}

