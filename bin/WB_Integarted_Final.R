```r
library(readr)
library(openxlsx)
library(stringr)
library(dplyr)

# ============================================================
# GET INPUT DIRECTORY FROM COMMAND LINE
#
# Usage:
# Rscript make_sv_workbook.R /path/to/input_directory
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: Rscript make_sv_workbook.R <input_directory>")
}

input_dir <- args[1]

if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
}

# ============================================================
# SAFE TSV READER
# ============================================================

read_tsv_safe <- function(file_path) {

    tryCatch(
        read_tsv(
            file_path,
            show_col_types = FALSE
        ),
        error = function(e) {
            warning(
                paste(
                    "Could not read:",
                    file_path
                )
            )
            return(NULL)
        }
    )
}

# ============================================================
# PICK ONE FILE FROM A LIST BASED ON PATTERN
# ============================================================

pick_file <- function(
    files,
    include_pattern,
    exclude_pattern = NULL
) {

    x <- files[
        grepl(
            include_pattern,
            basename(files),
            ignore.case = TRUE
        )
    ]

    if (!is.null(exclude_pattern)) {

        x <- x[
            !grepl(
                exclude_pattern,
                basename(x),
                ignore.case = TRUE
            )
        ]
    }

    x <- sort(x)

    if (length(x) == 0) {
        return(NULL)
    }

    return(x[1])
}

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <- file.path(
    input_dir,
    "Translocation_WB"
)

dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
)

message("Input directory: ", input_dir)
message("Output directory: ", output_dir)

# ============================================================
# CASES TO PROCESS
#
# NULL = process all case directories
#
# Example:
# cases_to_process <- c("CASE001", "CASE002")
# ============================================================

cases_to_process <- NULL

# ============================================================
# GET CASE DIRECTORIES
# ============================================================

case_dirs <- list.dirs(
    input_dir,
    recursive = FALSE,
    full.names = TRUE
)

# Don't treat output directory as a case
case_dirs <- case_dirs[
    basename(case_dirs) != "Translocation_WB"
]

if (!is.null(cases_to_process)) {

    case_dirs <- case_dirs[
        basename(case_dirs) %in% cases_to_process
    ]
}

# ============================================================
# PROCESS EACH CASE
# ============================================================

for (case_dir in case_dirs) {

    case <- basename(case_dir)

    message("")
    message("========================================")
    message("Processing case: ", case)
    message("========================================")

    # --------------------------------------------------------
    # Find all TSV files recursively
    # --------------------------------------------------------

    case_files <- list.files(
        case_dir,
        pattern = "\\.tsv$",
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(case_files) == 0) {

        warning(
            "No TSV files found for case: ",
            case
        )

        next
    }

    # --------------------------------------------------------
    # Find GRIDSS
    # --------------------------------------------------------

    gridss_file <- pick_file(
        case_files,
        include_pattern = "\\.gridss\\.tsv$",
        exclude_pattern = "somatic"
    )

    # --------------------------------------------------------
    # Find DELLY
    # --------------------------------------------------------

    delly_file <- pick_file(
        case_files,
        include_pattern = "_delly\\.tsv$",
        exclude_pattern = "_somatic"
    )

    # --------------------------------------------------------
    # Find SVABA
    # --------------------------------------------------------

    svaba_file <- pick_file(
        case_files,
        include_pattern = "svaba\\.somatic\\.sv\\.tsv$"
    )

    # --------------------------------------------------------
    # Find MANTA
    # --------------------------------------------------------

    manta_file <- pick_file(
        case_files,
        include_pattern = "manta_somaticSV\\.tsv$"
    )

    # --------------------------------------------------------
    # Find LUMPY
    # --------------------------------------------------------

    lumpy_file <- pick_file(
        case_files,
        include_pattern = "_svaba_lumpy\\.tsv$"
    )

    # --------------------------------------------------------
    # Read GRIDSS
    # --------------------------------------------------------

    gridss_df <- NULL

    if (!is.null(gridss_file)) {

        message("GRIDSS: ", gridss_file)

        gridss_df <- read_tsv_safe(
            gridss_file
        )
    }

    # --------------------------------------------------------
    # Read DELLY
    # --------------------------------------------------------

    delly_df <- NULL

    if (!is.null(delly_file)) {

        message("DELLY: ", delly_file)

        delly_df <- read_tsv_safe(
            delly_file
        )
    }

    # --------------------------------------------------------
    # Read SVABA
    # --------------------------------------------------------

    svaba_df <- NULL

    if (!is.null(svaba_file)) {

        message("SVABA: ", svaba_file)

        svaba_df <- read_tsv_safe(
            svaba_file
        )
    }

    # --------------------------------------------------------
    # Read MANTA
    # --------------------------------------------------------

    manta_df <- NULL

    if (!is.null(manta_file)) {

        message("MANTA: ", manta_file)

        manta_df <- read_tsv_safe(
            manta_file
        )
    }

    # --------------------------------------------------------
    # Read LUMPY
    # --------------------------------------------------------

    lumpy_df <- NULL

    if (!is.null(lumpy_file)) {

        message("LUMPY: ", lumpy_file)

        txt <- readLines(
            lumpy_file
        )

        # Find the line "lumpy output"
        start <- grep(
            "^lumpy output$",
            txt,
            ignore.case = TRUE
        )

        if (length(start) > 0) {

            # Use the first occurrence
            start <- start[1]

            # Read everything after "lumpy output"
            if (start < length(txt)) {

                lumpy_df <- read.delim(
                    textConnection(
                        txt[(start + 1):length(txt)]
                    ),
                    sep = "\t",
                    header = TRUE,
                    check.names = FALSE,
                    stringsAsFactors = FALSE
                )
            }

        } else {

            warning(
                "Could not find 'lumpy output' in: ",
                lumpy_file
            )
        }
    }

    # ========================================================
    # CREATE EXCEL WORKBOOK
    # ========================================================

    wb <- createWorkbook()

    if (!is.null(gridss_df)) {

        addWorksheet(
            wb,
            "GRIDSS"
        )

        writeData(
            wb,
            "GRIDSS",
            gridss_df
        )
    }

    if (!is.null(delly_df)) {

        addWorksheet(
            wb,
            "DELLY"
        )

        writeData(
            wb,
            "DELLY",
            delly_df
        )
    }

    if (!is.null(svaba_df)) {

        addWorksheet(
            wb,
            "SVABA"
        )

        writeData(
            wb,
            "SVABA",
            svaba_df
        )
    }

    if (!is.null(manta_df)) {

        addWorksheet(
            wb,
            "MANTA"
        )

        writeData(
            wb,
            "MANTA",
            manta_df
        )
    }

    if (!is.null(lumpy_df)) {

        addWorksheet(
            wb,
            "LUMPY"
        )

        writeData(
            wb,
            "LUMPY",
            lumpy_df
        )
    }

    # ========================================================
    # SAVE WORKBOOK
    # ========================================================

    if (length(wb$worksheets) > 0) {

        out_file <- file.path(
            output_dir,
            paste0(
                case,
                "_SV_calls.xlsx"
            )
        )

        saveWorkbook(
            wb,
            out_file,
            overwrite = TRUE
        )

        message(
            "Saved: ",
            out_file
        )

    } else {

        warning(
            "No matching SV caller TSVs found for case: ",
            case
        )
    }
}
```
