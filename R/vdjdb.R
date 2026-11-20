#' VDJdb data source
#'
#' VDJdb is a reference database of T-cell receptors with known antigen
#' specificities. `data_vdjdb()` manages a local copy of the slim VDJdb table
#' and can add matching VDJdb records to [immundata::ImmunData].
#'
#' @section Input and output:
#' `data <- data_vdjdb()` gives access to these commands:
#'
#' - **`data$download()`:** Print the download command and destination without
#'   downloading anything.
#' - **`data$download(execute = TRUE)`:** Download the latest VDJdb release.
#'   Use `force = TRUE` to replace an existing copy.
#' - **`data$check()`:** Report whether VDJdb has been downloaded.
#' - **`data$load()`:** Return the VDJdb table. The table stays on disk until
#'   you use `dplyr::collect()`. Use `data$load(output = "df")` to return a base
#'   R data frame instead.
#' - **`data$annotate(idata)`:** Match receptors in an `ImmunData` dataset to
#'   VDJdb and return a new dataset with VDJdb columns. By default, matching
#'   uses CDR3 amino-acid sequence and receptor locus.
#' - **`data$citation()`:** Return the VDJdb publication as an R citation.
#'
#' For the default matching, `ImmunData` must contain `junction_aa` and
#' `locus`. These are matched to the VDJdb columns `cdr3` and `gene`:
#' `c(junction_aa = "cdr3", locus = "gene")`. Use the `by` argument to choose
#' different matching columns. Added VDJdb columns start with `vdjdb_`, which
#' makes them easy to distinguish from the original columns.
#'
#' @param output_dir Directory used to store VDJdb. Defaults to the `data`
#'   subdirectory of `IMMUNDAPTER_DATA` when set, or of
#'   `tools::R_user_dir("immundapter", "cache")` otherwise.
#'
#' @return A VDJdb data source.
#'
#' @examplesIf identical(Sys.getenv("IMD_RUN_EXAMPLES"), "true")
#' data <- data_vdjdb(output_dir = tempfile("vdjdb-cache-"))
#'
#' # First call: print instructions. Second call: perform the download.
#' data$download()
#' data$download(execute = TRUE)
#' data$check()
#'
#' vdjdb_table <- data$load()
#' dplyr::collect(vdjdb_table)
#'
#' input_file <- system.file(
#'   "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
#'   package = "immundapter"
#' )
#' idata <- immundata::read_repertoires(
#'   path = input_file,
#'   schema = immundata::make_receptor_schema(
#'     features = c("v_call", "j_call", "junction_aa"), chains = "TRB"
#'   ),
#'   barcode_col = "cell_id", locus_col = "locus",
#'   umi_col = "duplicate_count",
#'   output_folder = tempfile("vdjdb-immundata-")
#' )
#' annotated_idata <- data$annotate(idata)
#' @seealso [VDJdb database](https://github.com/antigenomics/vdjdb-db)
#' @export
data_vdjdb <- function(output_dir = immundapter_data_path("data")) {
  VdjdbData$new(output_dir)
}


#' VDJdb data source
#'
#' @keywords internal
VdjdbData <- R6::R6Class(
  classname = "VdjdbData",
  inherit = AbstractData,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a VDJdb data source. Users should call [data_vdjdb()]
    #' instead.
    #' @param output_dir Directory used to store the VDJdb data file.
    initialize = function(output_dir = immundapter_data_path("data")) {
      super$initialize("vdjdb", output_dir)
    }
  ),
  private = list(
    download_impl = function(output_dir, execute, force) {
      path <- fs::path(output_dir, "vdjdb.slim.txt")

      if (!execute) {
        cli::cli_inform(c(
          "i" = "Download the latest VDJdb release:",
          ">" = "Run {.code data$download(execute = TRUE)}",
          "*" = "The database will be stored at {.path {path}}.",
          "*" = paste0(
            "Releases: {.url ",
            "https://github.com/antigenomics/vdjdb-db/releases}"
          )
        ))
        return(invisible(FALSE))
      }

      if (fs::file_exists(path) && !force) {
        cli::cli_inform(c(
          "v" = "VDJdb is already downloaded.",
          ">" = "Location: {.path {path}}",
          "i" = "Use {.code force = TRUE} to replace it."
        ))
        return(invisible(TRUE))
      }

      fs::dir_create(output_dir)
      manifest <- tempfile(fileext = ".txt")
      archive <- tempfile(fileext = ".zip")
      extract_dir <- tempfile()
      dir.create(extract_dir)
      on.exit(unlink(c(manifest, archive), force = TRUE), add = TRUE)
      on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)

      cli::cli_inform(c(
        "i" = "Downloading the latest VDJdb release..."
      ))
      utils::download.file(
        paste0(
          "https://raw.githubusercontent.com/antigenomics/",
          "vdjdb-db/master/latest-version.txt"
        ),
        manifest,
        quiet = TRUE,
        mode = "wb"
      )
      release_url <- readLines(manifest, n = 1L, warn = FALSE)
      utils::download.file(
        release_url,
        archive,
        quiet = TRUE,
        mode = "wb"
      )

      contents <- utils::unzip(archive, list = TRUE)$Name
      filenames <- basename(contents)
      entry <- contents[filenames == "vdjdb.slim.txt"]
      if (!length(entry)) {
        entry <- contents[filenames == "vdjdb.slim.txt.gz"]
      }
      if (length(entry) != 1L) {
        rlang::abort(
          "Expected one `vdjdb.slim.txt` or `vdjdb.slim.txt.gz` in the release."
        )
      }

      utils::unzip(archive, files = entry, exdir = extract_dir, junkpaths = TRUE)
      pending <- tempfile("vdjdb-", tmpdir = output_dir)
      on.exit(unlink(pending), add = TRUE)

      # gzfile() reads plain and gzip-compressed files.
      input <- gzfile(fs::path(extract_dir, basename(entry)), "rb")
      on.exit(close(input), add = TRUE, after = FALSE)
      output <- file(pending, "wb")
      tryCatch({
        repeat {
          bytes <- readBin(input, "raw", n = 1024L * 1024L)
          if (!length(bytes)) {
            break
          }
          writeBin(bytes, output)
        }
      }, finally = close(output))
      if (!file.rename(pending, path)) {
        rlang::abort("Failed to save the VDJdb database.")
      }

      cli::cli_inform(c(
        "v" = "VDJdb was downloaded.",
        ">" = "Location: {.path {path}}"
      ))
      invisible(TRUE)
    },
    check_impl = function(output_dir) {
      path <- fs::path(output_dir, "vdjdb.slim.txt")
      available <- isTRUE(fs::file_exists(path))

      if (available) {
        cli::cli_inform(c(
          "v" = "VDJdb is available.",
          ">" = "Location: {.path {path}}"
        ))
      } else {
        cli::cli_inform(c(
          "x" = "VDJdb is not downloaded.",
          ">" = "Expected location: {.path {path}}",
          "i" = "Run {.code data$download()} for instructions."
        ))
      }

      available
    },
    load_impl = function(output_dir) {
      path <- fs::path(output_dir, "vdjdb.slim.txt")
      if (!fs::file_exists(path)) {
        rlang::abort(
          "VDJdb is not downloaded. Run `data_vdjdb()$download()` first."
        )
      }

      duckplyr::read_csv_duckdb(
        path,
        options = list(
          delim = "\t",
          max_line_size = 10000000
        )
      )
    },
    annotate_impl = function(idata,
                             output_dir,
                             by = c(junction_aa = "cdr3", locus = "gene"),
                             keep_repertoires = TRUE) {
      checkmate::assert_character(by, min.len = 1L, names = "named")
      checkmate::assert_flag(keep_repertoires)

      database <- private$load_impl(output_dir)
      value_columns <- setdiff(colnames(database), unname(by))
      rename_columns <- rlang::set_names(
        rlang::syms(value_columns),
        paste0("vdjdb_", value_columns)
      )
      database <- database |>
        dplyr::rename(!!!rename_columns)

      immundata::annotate(
        idata = idata,
        annotations = database,
        by = by,
        keep_repertoires = keep_repertoires,
        remove_limit = TRUE
      )
    },
    citation_impl = function() {
      utils::bibentry(
        bibtype = "Article",
        title = paste(
          "VDJdb in the pandemic era: a compendium of T cell receptors",
          "specific for SARS-CoV-2"
        ),
        author = utils::as.person(c(
          "Mikhail Goncharov",
          "Dmitry Bagaev",
          "Dmitrii Shcherbinin",
          "Ivan Zvyagin",
          "Dmitry Bolotin",
          "Paul G. Thomas",
          "Anastasia A. Minervina",
          "Mikhail V. Pogorelyy",
          "Kristin Ladell",
          "James E. McLaren",
          "David A. Price",
          "Thi H. O. Nguyen",
          "Louise C. Rowntree",
          "E. Bridie Clemens",
          "Katherine Kedzierska",
          "Garry Dolton",
          "Cristina Rafael Rius",
          "Andrew Sewell",
          "Jerome Samir",
          "Fabio Luciani",
          "Ksenia V. Zornikova",
          "Alexandra A. Khmelevskaya",
          "Saveliy A. Sheetikov",
          "Grigory A. Efimov",
          "Dmitry Chudakov",
          "Mikhail Shugay"
        )),
        journal = "Nature Methods",
        year = "2022",
        volume = "19",
        pages = "1017--1019",
        doi = "10.1038/s41592-022-01578-0"
      )
    }
  )
)
