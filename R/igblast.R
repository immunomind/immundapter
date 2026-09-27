#' @title Assign germline genes to raw receptor sequences with IgBLAST
#'
#' @description
#'
#' IgBLAST identifies the germline genes that most likely produced each immune
#' receptor sequence. It produces V, D, J, and, when available, C gene calls
#' in an AIRR table that can be imported with [immundata::read_repertoires()].
#'
#' The usual workflow is:
#'
#' ```r
#' run <- tool_igblast()$run(c(
#'   fasta = "sequences.fasta",
#'   annotations = "contig_annotations.csv"
#' ))
#' idata <- immundata::read_repertoires(
#'   path = run$results()$final_output,
#'   schema = immundata::make_receptor_schema(
#'     features = c("v_call", "j_call", "junction_aa"),
#'     chains = c("IGH", "IGK|IGL")
#'   ),
#'   barcode_col = "barcode",
#'   locus_col = "locus",
#'   umi_col = "umis",
#'   preprocess = immundata::make_default_preprocessing("10x"),
#'   rename_columns = NULL
#' )
#' ```
#'
#' @details
#' IgBLAST compares observed receptor sequences with a germline database. The
#' database selected through `igblastr` determines the species and receptor
#' type. For example, select a human BCR database before analysing human
#' immunoglobulin sequences.
#'
#' Query sequences must be **nucleotide sequences without alignment gaps**.
#' Supply the complete assembled V(D)J sequence when possible. Amino-acid
#' sequences and sequences that already contain IMGT dots are not valid input.
#' A CDR3 nucleotide sequence alone usually does not contain enough information
#' for useful V- and J-gene annotation.
#'
#' @section Before the first analysis:
#'
#' Install IgBLAST and one germline database with `igblastr`. Select the
#' database in the same R session in which you run the analysis:
#'
#' ```r
#' igblastr::install_igblast()
#' igblastr::install_IMGT_germline_db(
#'   "<release>",
#'   organism = "Homo sapiens",
#'   tcr.db = FALSE
#' )
#' igblastr::use_germline_db("<db-name>")
#' ```
#'
#' Use `tool$help()` for setup commands and `tool$check()` to confirm that the
#' selected database and IgBLAST installation are ready.
#'
#' @section Run the tool:
#'
#' `tool_igblast()` creates a reusable tool. Its main methods are:
#'
#' * `tool$help()`: print installation and database setup instructions;
#' * `tool$check()`: check that IgBLAST and a germline database are ready;
#' * `tool$run(input, config = NULL, output_dir = NULL)`: run one analysis.
#'
#' The tool keeps its default settings in `tool$config`. Supplying `config` to
#' `$run()` replaces those settings for that run only.
#'
#' @section Choose the input:
#'
#' `tool$run()` accepts two forms of input:
#'
#' * **A FASTA file.** Pass a path ending in `.fasta` or `.fasta.gz`. Use this
#'   when no annotation table needs to be preserved.
#' * **A FASTA file and its annotation table.** Pass
#'   `c(fasta = "sequences.fasta", annotations = "annotations.csv")`. Use this
#'   for data such as 10x contigs when barcodes, UMI counts, and sample columns
#'   must remain beside the new IgBLAST annotations.
#'
#' @section Read the results:
#'
#' `tool$run()` saves the analysis and returns a completed run object. Its
#' `$path` is the output directory and `$artifacts` lists the saved files. Call
#' `run$results()` to access the biological results. Every run provides:
#'
#' * `results$airr`: the complete AIRR annotation table;
#' * `results$mapping`: the link between input records and AIRR `sequence_id`;
#' * `results$raw_output`: the path to the saved AIRR TSV file.
#'
#' The tables are disk-backed. Use `dplyr::collect()` only when you want to
#' bring a table into R memory.
#'
#' With paired FASTA and annotation input, `results$final` contains the source
#' columns joined with all AIRR columns. `results$final_output` is the path to
#' the same table as a TSV file. Pass this file to
#' [immundata::read_repertoires()] to create an `ImmunData` dataset. Choose
#' preprocessing that retains `sequence_alignment` when subsequent analyses
#' need IMGT-gapped alignments and coordinates.
#'
#' @section IMGT gaps and coordinates:
#'
#' `igblastr` builds the searchable IgBLAST database from **ungapped** germline
#' sequences. It also keeps the original **IMGT-gapped V references**. These
#' two forms have different V-gene coordinates.
#'
#' For a selected IMGT database, immundapter automatically uses both forms from
#' the **same installed database release**. IgBLAST first performs the search
#' against the ungapped database. immundapter then uses the matching gapped V
#' allele to place IMGT dots in `sequence_alignment` and convert
#' `v_germline_start` and `v_germline_end` to the IMGT-gapped coordinate system.
#' This follows the normalization used by Change-O for germline reconstruction.
#'
#' Do not download a newer gapped FASTA separately: allele content and gap
#' positions must match the database used by IgBLAST. If the selected IMGT
#' database does not contain its original gapped V FASTA files, reinstall that
#' database with `igblastr::install_IMGT_germline_db()`.
#'
#' Normalization uses the first allele in `v_call`. If the allele or the exact
#' aligned V end cannot be matched to the selected reference, that row keeps
#' the original IgBLAST alignment and coordinates. Non-IMGT databases are also
#' returned in their original IgBLAST coordinate system. D- and J-region
#' coordinates are not converted to IMGT-gapped V coordinates.
#'
#' @section Input details:
#'
#' FASTA names must be unique. In `>cell_1 description`, the name is `cell_1`.
#' For paired input, the annotation table must contain a unique `contig_id` or
#' `sequence_id` column whose values match the FASTA names. CSV, TSV, TXT, and
#' gzip-compressed annotation files are accepted.
#'
#' @param num_threads Number of computer threads available to IgBLAST. `NULL`
#'   uses the IgBLAST default.
#' @param run_args Additional named arguments passed through `...` in
#'   `igblastr::igblastn()`. This is intended for advanced settings. Arguments
#'   required by the immundapter workflow cannot be changed here.
#' @param config Settings created by [conf_igblast()].
#'
#' @return [conf_igblast()] returns read-only IgBLAST settings.
#'   [tool_igblast()] returns a reusable IgBLAST tool.
#'
#' @examplesIf identical(Sys.getenv("IMD_RUN_EXAMPLES"), "true")
#' # Select an installed human BCR IMGT database.
#' databases <- igblastr::list_germline_dbs(names.only = TRUE)
#' databases <- grep(
#'   "Homo_sapiens.IGH+IGK+IGL",
#'   databases,
#'   value = TRUE,
#'   fixed = TRUE
#' )
#' databases <- grep("^IMGT-", databases, value = TRUE)
#' stopifnot(length(databases) > 0L)
#' suppressMessages(
#'   igblastr::use_germline_db(sort(databases)[[length(databases)]])
#' )
#'
#' config <- conf_igblast()
#' tool <- tool_igblast(config)
#' stopifnot(tool$check())
#'
#' antitnf_dir <- system.file("extdata", "antitnf", package = "immundapter")
#' antitnf_fasta <- file.path(antitnf_dir, "bcr_contigs.fasta.gz")
#' antitnf_annotations <- file.path(
#'   antitnf_dir,
#'   "contig_annotations.csv.gz"
#' )
#'
#' # Analyse a FASTA file directly.
#' fasta_run <- tool$run(
#'   antitnf_fasta,
#'   output_dir = tempfile("igblast-fasta-")
#' )
#' fasta_airr <- dplyr::collect(fasta_run$results()$airr)
#'
#' # Preserve the matching 10x annotation table and import the AIRR result.
#' paired_run <- tool$run(
#'   c(
#'     fasta = antitnf_fasta,
#'     annotations = antitnf_annotations
#'   ),
#'   output_dir = tempfile("igblast-10x-")
#' )
#' idata <- immundata::read_repertoires(
#'   path = paired_run$results()$final_output,
#'   schema = immundata::make_receptor_schema(
#'     features = c("v_call", "j_call", "junction_aa"),
#'     chains = c("IGH", "IGK|IGL")
#'   ),
#'   barcode_col = "barcode",
#'   locus_col = "locus",
#'   umi_col = "umis",
#'   preprocess = immundata::make_default_preprocessing("10x"),
#'   rename_columns = NULL,
#'   output_folder = tempfile("antitnf-immundata-")
#' )
#' @name tool_igblast
#' @seealso
#'   [igblastr overview](https://bioconductor.org/packages/release/bioc/vignettes/igblastr/inst/doc/igblastr_overview.html)
#' @export
conf_igblast <- function(num_threads = NULL,
                         run_args = list()) {
  IgblastToolConfig$new(list(
    num_threads = num_threads,
    run_args = run_args
  ))
}


#' @rdname tool_igblast
#' @export
tool_igblast <- function(config = conf_igblast()) {
  IgblastTool$new(config)
}


#' @noRd
IgblastToolConfig <- R6::R6Class(
  classname = "IgblastToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    num_threads = NULL,
    run_args = NULL,
    #' @description
    #' Create checked IgBLAST settings. Users should call
    #' [conf_igblast()] instead.
    #' @param values A named list of IgBLAST settings.
    initialize = function(values) {
      checkmate::assert_int(values$num_threads, lower = 1L, null.ok = TRUE)
      checkmate::assert_list(values$run_args)

      argument_names <- names(values$run_args)
      if (
        length(values$run_args) > 0L &&
          (
            is.null(argument_names) ||
              any(!nzchar(argument_names)) ||
              anyDuplicated(argument_names)
          )
      ) {
        rlang::abort("IgBLAST `run_args` must have unique, non-empty names.")
      }

      fixed <- c(
        "query", "outfmt", "germline_db_V", "germline_db_D",
        "germline_db_J", "organism", "c_region_db",
        "custom_internal_data", "auxiliary_data", "domain_system",
        "ig_seqtype", "out", "parse.out", "show.in.browser",
        "show.command.only", "num_threads"
      )
      overridden <- intersect(argument_names, fixed)
      if (length(overridden) > 0L) {
        rlang::abort(sprintf(
          "IgBLAST `run_args` cannot override: %s.",
          paste(overridden, collapse = ", ")
        ))
      }

      super$initialize(values)
    }
  )
)


#' @noRd
IgblastTool <- R6::R6Class(
  classname = "IgblastTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an IgBLAST tool. Users should call [tool_igblast()] instead.
    #' @param config Settings returned by [conf_igblast()].
    initialize = function(config = conf_igblast()) {
      checkmate::assert_r6(config, classes = "IgblastToolConfig")
      super$initialize("igblast", config, IgblastRun)
    }
  ),
  private = list(
    help_impl = function() {
      cli::cli_inform(c(
        "i" = "Set up IgBLAST through {.pkg igblastr}:",
        ">" = '{.code pak::pak("bioc::igblastr")}',
        ">" = "{.code igblastr::install_igblast()}",
        ">" = "{.code igblastr::list_germline_dbs()}",
        ">" = paste(
          "{.code igblastr::install_IMGT_germline_db(",
          '"<release>", organism = "Homo sapiens", tcr.db = TRUE)}'
        ),
        ">" = '{.code igblastr::use_germline_db("<db-name>")}',
        "i" = paste(
          "The selected database determines organism and Ig versus TCR.",
          "Selection must happen in the R session that calls `$run()`."
        ),
        "i" = paste(
          "A compatible bundled IMGT C-region database is selected",
          "automatically when no C-region database is already selected."
        )
      ))
    },
    check_impl = function() {
      runtime <- igblast_runtime_info()
      if (!runtime$package_installed) {
        cli::cli_inform(c(
          "x" = "{.pkg igblastr} is not installed.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
        return(FALSE)
      }
      if (!runtime$igblast_installed) {
        cli::cli_inform(c(
          "v" = "{.pkg igblastr} {runtime$igblastr_version} is installed.",
          "x" = "No working IgBLAST installation was found."
        ))
        return(FALSE)
      }
      if (is.null(runtime$germline_db)) {
        cli::cli_inform(c(
          "v" = "IgBLAST {runtime$igblast_version} is available.",
          "x" = "No germline V/D/J database is selected in this R session."
        ))
        return(FALSE)
      }

      automatic_c_region_db <- if (nzchar(runtime$c_region_db)) {
        ""
      } else {
        igblast_find_bundled_c_region_db(runtime$germline_db)
      }

      cli::cli_inform(c(
        "v" = "{.pkg igblastr} {runtime$igblastr_version} is installed.",
        "v" = "IgBLAST {runtime$igblast_version} is available.",
        "v" = "Selected germline database: {runtime$germline_db}",
        "i" = if (nzchar(runtime$c_region_db)) {
          "Selected C-region database: {runtime$c_region_db}"
        } else if (nzchar(automatic_c_region_db)) {
          paste(
            "Bundled C-region database selected automatically at run time:",
            "{automatic_c_region_db}"
          )
        } else {
          paste(
            "No selected or compatible bundled C-region database was found;",
            "IgBLAST will not assign `c_call`."
          )
        }
      ))
      TRUE
    },
    citation_impl = function() {
      rlang::abort("IgBLAST `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      checkmate::assert_character(
        input,
        min.len = 1L,
        max.len = 2L,
        any.missing = FALSE
      )
      if (length(input) == 2L) {
        if (
          is.null(names(input)) ||
            !setequal(names(input), c("fasta", "annotations")) ||
            anyDuplicated(names(input))
        ) {
          rlang::abort(paste(
            "Paired IgBLAST input must be a character vector named",
            "`fasta` and `annotations`."
          ))
        }
        input <- input[c("fasta", "annotations")]
      } else {
        checkmate::assert_string(input, min.chars = 1L)
      }
      checkmate::assert_file_exists(input)
      input_names <- names(input)
      input <- fs::path_abs(input)
      names(input) <- input_names
      input
    },
    run_impl = function(input, config, run_dir) {
      runtime <- igblast_runtime_info()
      if (!runtime$package_installed || !runtime$igblast_installed) {
        rlang::abort(
          "igblastr and a working IgBLAST installation are required."
        )
      }
      if (is.null(runtime$germline_db)) {
        rlang::abort(paste(
          "No germline V/D/J database is selected.",
          "Call `igblastr::use_germline_db()` in this R session."
        ))
      }

      if (!nzchar(runtime$c_region_db)) {
        c_region_db <- igblast_find_bundled_c_region_db(runtime$germline_db)
        if (nzchar(c_region_db)) {
          selected <- tryCatch(
            {
              suppressMessages(
                getExportedValue("igblastr", "use_c_region_db")(c_region_db)
              )
              TRUE
            },
            error = function(error) {
              cli::cli_warn(c(
                "Unable to prepare the bundled C-region database",
                "{.val {c_region_db}}: {conditionMessage(error)}"
              ))
              FALSE
            }
          )
          if (selected) {
            on.exit(
              suppressMessages(
                getExportedValue("igblastr", "use_c_region_db")("")
              ),
              add = TRUE
            )
            runtime <- igblast_runtime_info()
          }
        }
      }

      prepared <- igblast_prepare_input(input)
      input_path <- fs::path(run_dir, "igblast-input.fasta")
      mapping_path <- fs::path(run_dir, "igblast-query-map.parquet")
      output_path <- fs::path(run_dir, "igblast-airr.tsv")
      raw_output_path <- tempfile("igblast-airr-", fileext = ".tsv")
      on.exit(unlink(raw_output_path), add = TRUE)
      final_path <- fs::path(run_dir, "igblast-final.tsv")
      config_path <- fs::path(run_dir, "igblast-config.json")

      writeLines(
        as.vector(rbind(
          paste0(">", prepared$queries$igblast_query_id),
          prepared$queries$sequence
        )),
        input_path,
        useBytes = TRUE
      )
      duckplyr::compute_parquet(prepared$mapping, mapping_path)

      fixed <- list(
        query = as.character(fs::path_abs(input_path)),
        outfmt = "AIRR",
        germline_db_V = "auto",
        germline_db_D = "auto",
        germline_db_J = "auto",
        organism = "auto",
        c_region_db = "auto",
        custom_internal_data = "auto",
        auxiliary_data = "auto",
        domain_system = "imgt",
        ig_seqtype = "auto",
        out = as.character(fs::path_abs(raw_output_path)),
        parse.out = FALSE,
        show.in.browser = FALSE,
        show.command.only = FALSE
      )
      arguments <- c(
        fixed,
        if (is.null(config$num_threads)) {
          list()
        } else {
          list(num_threads = config$num_threads)
        },
        config$run_args
      )

      igblast_call(arguments)
      if (!fs::file_exists(raw_output_path)) {
        rlang::abort("IgBLAST did not produce AIRR output.")
      }

      airr <- duckplyr::read_csv_duckdb(
        raw_output_path,
        options = list(
          delim = "\t",
          header = TRUE,
          types = list(c(sequence_id = "VARCHAR"))
        )
      )
      igblast_validate_airr(airr, prepared$queries)
      v_references <- igblast_find_imgt_v_references(runtime$germline_db)
      if (length(v_references)) {
        airr <- igblast_normalize_imgt(
          airr,
          igblast_make_imgt_v_plans(v_references)
        )
      }
      duckplyr::compute_csv(
        airr,
        output_path,
        options = list(delimiter = "\t", header = TRUE)
      )

      artifacts <- list(
        input = input_path,
        config = config_path,
        mapping = mapping_path,
        airr = output_path
      )
      if (!is.null(prepared$annotations)) {
        final <- igblast_merge_annotations(
          airr,
          prepared$annotations,
          prepared$annotation_id_col
        )
        duckplyr::compute_csv(
          final,
          final_path,
          options = list(delimiter = "\t", header = TRUE)
        )
        artifacts$final <- final_path
      }

      input_count <- prepared$mapping |>
        dplyr::summarise(number_input_sequences = dplyr::n()) |>
        dplyr::collect() |>
        as.data.frame()

      jsonlite::write_json(
        list(
          input_type = if (length(input) == 2L) {
            "fasta_annotations"
          } else {
            "fasta"
          },
          num_threads = config$num_threads,
          fixed = fixed[setdiff(names(fixed), c("query", "out"))],
          run_args = config$run_args,
          effective_args = arguments[
            setdiff(names(arguments), c("query", "out"))
          ],
          runtime = runtime,
          number_input_sequences = input_count$number_input_sequences,
          number_unique_sequences = nrow(prepared$queries)
        ),
        path = config_path,
        auto_unbox = TRUE,
        null = "null",
        pretty = TRUE
      )

      artifacts
    }
  )
)
