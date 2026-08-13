#' Assign immune-receptor germline genes with IgBLAST
#'
#' IgBLAST compares immune-receptor nucleotide sequences with germline gene
#' databases. It reports V, D, J, and optional constant-region assignments in
#' AIRR format. The database selected with `igblastr::use_germline_db()`
#' determines the organism and whether the input contains B-cell or T-cell
#' receptors.
#'
#' @section Input to `tool$run()`:
#' Choose one of these three input forms:
#'
#' - **`ImmunData`:** Pass an [immundata::ImmunData] dataset. `sequence_col` is
#'   the name of the column with full-length, ungapped nucleotide sequences.
#'   Missing and empty sequences are skipped. Exact duplicate sequences are
#'   analyzed once and then matched back to all original chains.
#' - **FASTA file:** Pass one file path, for example `"sequences.fasta"` or
#'   `"sequences.fasta.gz"`. Every FASTA sequence must have a unique name.
#'   Each sequence must contain ungapped IUPAC nucleotide letters. Every record is
#'   analyzed separately, even when two records have the same sequence.
#' - **FASTA file with an annotation table:** Pass two named paths:
#'   `c(fasta = "sequences.fasta", annotations = "annotations.csv")`. The
#'   annotation table may be CSV, TSV, or tab-delimited TXT, with optional gzip
#'   compression. It usually comes from 10x Genomics. It must contain a unique
#'   `contig_id` or `sequence_id` column. Values in this column must match the
#'   FASTA sequence names exactly. This form keeps barcodes, read counts, UMI
#'   counts, and other source columns while re-annotating 10x Genomics data.
#'
#' In a FASTA header such as `>cell_1 description`, the sequence name is
#' `cell_1`. IgBLAST writes this name to the AIRR `sequence_id` column.
#'
#' @section Output from `tool$run()`:
#' `tool$run()` returns an [IgblastRun]. It records where the results were
#' saved. To read the result tables, first run `results <- run$results()`.
#'
#' [IgblastResults] provides:
#'
#' - **`results$airr`:** The complete AIRR table. It stays on disk until you
#'   use `dplyr::collect(results$airr)`.
#' - **`results$mapping`:** A table that shows which input chain or FASTA
#'   sequence produced each AIRR `sequence_id`.
#' - **`results$raw_output`:** The file path of the AIRR TSV result.
#'
#' Additional output depends on the input form:
#'
#' - **`ImmunData`:** `results$annotate(idata)` returns a new `ImmunData`
#'   dataset with selected AIRR fields added to its chain table. By
#'   default, it adds `sequence_alignment` and `germline_alignment`.
#' - **FASTA file:** Use `$airr`, `$mapping`, and `$raw_output`. These results
#'   cannot add columns to `ImmunData` because the input did not come from an
#'   `ImmunData` dataset.
#' - **FASTA file with an annotation table:** `results$final` contains the
#'   original annotation table joined with the AIRR columns.
#'   `results$final_output` is the path to the same table as a TSV file, ready
#'   to read into `ImmunData`.
#'
#' @param sequence_col For `ImmunData` input, the name of the chain column
#'   containing full-length, ungapped nucleotide sequences. This setting is not
#'   used for FASTA input.
#' @param num_threads Number of computer threads used by IgBLAST. `NULL` uses
#'   the IgBLAST default.
#' @param run_args Additional named arguments passed through `...` in
#'   `igblastr::igblastn()`. Use this only for advanced IgBLAST settings.
#'   Arguments controlled by immundapter cannot be changed here.
#' @param config An [IgblastToolConfig] created by [conf_igblast()].
#'
#' @return [conf_igblast()] returns a read-only [IgblastToolConfig].
#'   [tool_igblast()] returns an IgBLAST tool.
#'
#' @examplesIf interactive()
#' config <- conf_igblast(sequence_col = "sequence")
#' tool <- tool_igblast(config)
#' tool$help()
#'
#' # Installation
#' pak::pak("bioc::igblastr")
#' igblastr::install_igblast()
#' igblastr::install_IMGT_germline_db(
#'   "<release>",
#'   organism = "Homo sapiens",
#'   tcr.db = TRUE
#' )
#' igblastr::use_germline_db("<db-name>")
#'
#' # Workflow 1: add IgBLAST results to ImmunData.
#' run <- tool$run(idata, output_dir = "analysis/igblast")
#' results <- run$results()
#' dplyr::collect(results$airr)
#' dplyr::collect(results$mapping)
#' aligned_idata <- results$annotate(idata)
#'
#' # Workflow 2: analyze a FASTA file.
#' fasta_run <- tool$run("queries.fasta", output_dir = "analysis/igblast")
#' fasta_airr <- dplyr::collect(fasta_run$results()$airr)
#'
#' # Workflow 3: analyze FASTA and keep a matching annotation table.
#' paired_run <- tool$run(
#'   c(
#'     fasta = "filtered_contig.fasta",
#'     annotations = "filtered_contig_annotations.csv"
#'   ),
#'   output_dir = "analysis/igblast-10x"
#' )
#' paired_results <- paired_run$results()
#' dplyr::collect(paired_results$final)
#' paired_results$final_output
#'
#' bcr_schema <- immundata::make_receptor_schema(
#'   features = c("v_call", "j_call", "junction_aa"),
#'   chains = c("IGH", "IGK|IGL")
#' )
#' bcr_data <- immundata::read_repertoires(
#'   path = paired_results$final_output,
#'   schema = bcr_schema,
#'   barcode_col = "barcode",
#'   locus_col = "locus",
#'   umi_col = "umis",
#'   preprocess = immundata::make_default_preprocessing("10x"),
#'   postprocess = NULL,
#'   rename_columns = NULL,
#'   output_folder = "analysis/immundata-bcr"
#' )
#' @name igblast
#' @seealso
#'   [igblastr overview](https://bioconductor.org/packages/release/bioc/vignettes/igblastr/inst/doc/igblastr_overview.html)
#' @export
conf_igblast <- function(sequence_col = "sequence",
                         num_threads = NULL,
                         run_args = list()) {
  IgblastToolConfig$new(list(
    sequence_col = sequence_col,
    num_threads = num_threads,
    run_args = run_args
  ))
}


#' @rdname igblast
#' @export
tool_igblast <- function(config = conf_igblast()) {
  IgblastTool$new(config)
}


#' @rdname igblast
#'
#' @description
#' `IgblastToolConfig` stores validated IgBLAST settings. It contains no
#' sequences and no results. Create it with [conf_igblast()] and pass it to
#' [tool_igblast()]. The settings cannot be changed after creation. Call
#' [conf_igblast()] again when different settings are needed.
#'
#' @field sequence_col Column containing ungapped nucleotide query sequences.
#' @field num_threads Optional number of IgBLAST threads.
#' @field run_args Additional named IgBLAST command-line arguments.
IgblastToolConfig <- R6::R6Class(
  classname = "IgblastToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    sequence_col = NULL,
    num_threads = NULL,
    run_args = NULL,
    #' @description
    #' Create checked IgBLAST settings. Users should call
    #' [conf_igblast()] instead.
    #' @param values A named list of IgBLAST settings.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_string(values$sequence_col)
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


#' Run IgBLAST analyses
#'
#' `IgblastTool` uses one [IgblastToolConfig] as its default settings. Create it
#' with [tool_igblast()]. See [igblast] for accepted input and available output.
#'
#' @section Methods:
#'
#' - **`$help()`:** Print IgBLAST and germline-database setup instructions.
#' - **`$check()`:** Report whether igblastr, IgBLAST, and a germline database
#'   are ready. Return `TRUE` when all required components are available.
#' - **`$run(input, config = NULL, output_dir = NULL)`:** Run IgBLAST and return
#'   an [IgblastRun]. See [igblast] for the three accepted input forms.
#' - **`$citation()`:** Reserved for future citation support.
#'
#' @keywords internal
IgblastTool <- R6::R6Class(
  classname = "IgblastTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an IgBLAST tool. Users should call [tool_igblast()] instead.
    #' @param config An [IgblastToolConfig] returned by [conf_igblast()].
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
          "A C-region database is optional. Install and select one with",
          "{.code install_IMGT_c_region_db()} and",
          "{.code use_c_region_db()} when constant-region calls are needed."
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

      cli::cli_inform(c(
        "v" = "{.pkg igblastr} {runtime$igblastr_version} is installed.",
        "v" = "IgBLAST {runtime$igblast_version} is available.",
        "v" = "Selected germline database: {runtime$germline_db}",
        "i" = if (nzchar(runtime$c_region_db)) {
          "Selected C-region database: {runtime$c_region_db}"
        } else {
          "No C-region database is selected."
        }
      ))
      TRUE
    },
    citation_impl = function() {
      rlang::abort("IgBLAST `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      if (inherits(input, "ImmunData")) {
        checkmate::assert_r6(input, classes = "ImmunData")
        return(input)
      }
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

      prepared <- igblast_prepare_input(input, config$sequence_col)
      input_path <- fs::path(run_dir, "igblast-input.fasta")
      mapping_path <- fs::path(run_dir, "igblast-query-map.parquet")
      output_path <- fs::path(run_dir, "igblast-airr.tsv")
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
        out = as.character(fs::path_abs(output_path)),
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

      jsonlite::write_json(
        list(
          input_type = if (inherits(input, "ImmunData")) {
            "immundata"
          } else if (length(input) == 2L) {
            "fasta_annotations"
          } else {
            "fasta"
          },
          sequence_col = config$sequence_col,
          num_threads = config$num_threads,
          run_args = config$run_args,
          fixed = fixed[setdiff(names(fixed), c("query", "out"))],
          runtime = runtime
        ),
        path = config_path,
        auto_unbox = TRUE,
        null = "null",
        pretty = TRUE
      )

      igblast_call(arguments)
      if (!fs::file_exists(output_path)) {
        rlang::abort("IgBLAST did not produce AIRR output.")
      }

      artifacts <- list(
        input = input_path,
        config = config_path,
        mapping = mapping_path,
        airr = output_path
      )
      if (!is.null(prepared$annotations)) {
        airr <- duckplyr::read_csv_duckdb(
          output_path,
          options = list(delim = "\t", header = TRUE)
        )
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

      artifacts
    }
  )
)
