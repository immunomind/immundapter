#' VDJdb data source
#'
#' `data_vdjdb()` provides access to the slim VDJdb table for single-chain TCR
#' annotation.
#'
#' Use `$download()` to print download instructions, or
#' `$download(execute = TRUE)` to download the latest VDJdb release. `$check()`
#' reports whether the data are available and where they are stored. `$load()`
#' returns the table as a duckplyr data frame.
#'
#' `$annotate()` joins VDJdb to an [immundata::ImmunData] object by CDR3 amino
#' acid sequence and locus by default:
#' `c(junction_aa = "cdr3", locus = "gene")`. Non-key VDJdb columns are prefixed
#' with `vdjdb_`.
#'
#' @param output_dir Directory used to store VDJdb. Defaults to the `data`
#'   directory under `IMMUNDAPTER_DATA`, or `~/immundapter_data/data` when the
#'   environment variable is not set.
#'
#' @return A VDJdb data-source adapter.
#'
#' @seealso \url{https://github.com/antigenomics/vdjdb-db}
#' @export
data_vdjdb <- function(output_dir = immundapter_data_path("data")) {
  VdjdbData$new(output_dir)
}


#' VDJdb data-source adapter
#'
#' @keywords internal
VdjdbData <- R6::R6Class(
  classname = "VdjdbData",
  inherit = AbstractData,
  cloneable = FALSE,
  public = list(
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
      entry <- contents[endsWith(contents, "/vdjdb.slim.txt")]
      if (length(entry) != 1L) {
        rlang::abort("The VDJdb release does not contain `vdjdb.slim.txt`.")
      }

      utils::unzip(archive, files = entry, exdir = extract_dir)
      if (!file.copy(
        fs::path(extract_dir, entry),
        path,
        overwrite = TRUE
      )) {
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
