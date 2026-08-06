library(readr)
library(openxlsx)
library(stringr)
library(dplyr)

# SAFE TSV READER

read_tsv_safe <- function(file_path) {
  tryCatch(
    read_tsv(file_path, show_col_types = FALSE),
    error = function(e) {
      warning(paste("Could not read:", file_path))
      return(NULL)
    }
  )
}

# PICK ONE FILE FROM A LIST BASED ON PATTERN

pick_file <- function(files, include_pattern, exclude_pattern = NULL) {
  x <- files[grepl(include_pattern, basename(files), ignore.case = TRUE)]
  
  if (!is.null(exclude_pattern)) {
    x <- x[!grepl(exclude_pattern, basename(x), ignore.case = TRUE)]
  }
  
  x <- sort(x)
  
  if (length(x) == 0) return(NULL)
  return(x[1])
}

# Input/Output

input_dir <- ""

output_dir <- file.path(input_dir, "Translocation_WB")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cases_to_process <- NULL

# Get case folders

case_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
case_dirs <- case_dirs[basename(case_dirs) != "Translocation_WB"]

if (!is.null(cases_to_process)) {
  case_dirs <- case_dirs[basename(case_dirs) %in% cases_to_process]
}

# Process files with .tsv extension

for (case_dir in case_dirs) {
  
  case <- basename(case_dir)
  message("Processing case: ", case)
  
  # Get all TSV files inside this case folder and subfolders
  case_files <- list.files(
    case_dir,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(case_files) == 0) {
    warning("No TSV files found for case: ", case)
    next
  }
  
  
  # Annotation of the file required 
  
  gridss_file <- pick_file(
    case_files,
    include_pattern = "\\.gridss\\.tsv$",
    exclude_pattern = "somatic"
  )
  
  delly_file <- pick_file(
    case_files,
    include_pattern = "_delly\\.tsv$",
    exclude_pattern = "_somatic"
  )
  
  svaba_file <- pick_file(
    case_files,
    include_pattern = "svaba\\.somatic\\.sv\\.tsv$"
  )
  
  manta_file <- pick_file(
    case_files,
    include_pattern = "manta_somaticSV\\.tsv$"
  )
  
  lumpy_file <- pick_file(
    case_files,
    include_pattern = "_svaba_lumpy\\.tsv$"
  )
  
  
  # Read selected files
  
  
  gridss_df <- if (!is.null(gridss_file)) read_tsv_safe(gridss_file) else NULL
  delly_df  <- if (!is.null(delly_file))  read_tsv_safe(delly_file)  else NULL
  svaba_df  <- if (!is.null(svaba_file))  read_tsv_safe(svaba_file)  else NULL
  manta_df  <- if (!is.null(manta_file))  read_tsv_safe(manta_file)  else NULL
  lumpy_df <- NULL
  
  if (!is.null(lumpy_file)) {
    
    txt <- readLines(lumpy_file)
    
    # Find "lumpy output"
    start <- grep("^lumpy output$", txt, ignore.case = TRUE)
    
    if (length(start) > 0) {
      
      # Read only from the LUMPY header onwards
      lumpy_df <- read.delim(
        textConnection(txt[(start + 1):length(txt)]),
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
    }
    
  }
  
  # Create workbook
  
  wb <- createWorkbook()
  
  if (!is.null(gridss_df)) {
    addWorksheet(wb, "GRIDSS")
    writeData(wb, "GRIDSS", gridss_df)
  }
  
  if (!is.null(delly_df)) {
    addWorksheet(wb, "DELLY")
    writeData(wb, "DELLY", delly_df)
  }
  
  if (!is.null(svaba_df)) {
    addWorksheet(wb, "SVABA")
    writeData(wb, "SVABA", svaba_df)
  }
  
  if (!is.null(manta_df)) {
    addWorksheet(wb, "MANTA")
    writeData(wb, "MANTA", manta_df)
  }
  
  if (!is.null(lumpy_df)) {
    addWorksheet(wb, "LUMPY")
    writeData(wb, "LUMPY", lumpy_df)
  }
  
  # Save workbook if at least one caller exists
  
   if (length(wb$worksheets) > 0) {
    out_file <- file.path(output_dir, paste0(case, "_SV_calls.xlsx"))
    saveWorkbook(wb, out_file, overwrite = TRUE)
    message("Saved: ", out_file)
  } else {
    warning("No matching SV caller TSVs found for case: ", case)
  }
}

