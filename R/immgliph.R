#' Group similar T-cell receptors with ImmGLIPH
#'
#' ImmGLIPH groups T-cell receptor beta chains that share sequence patterns and
#' may recognize related antigens. The analysis uses CDR3 amino-acid sequences
#' and can also use V genes, donor names, HLA types, and clone counts.
#'
#' @section Input to `tool$run()`:
#' `tool$run()` accepts one [immundata::ImmunData] dataset. It does not accept a
#' text file directly. The dataset must contain TRB receptors. For paired
#' TRA/TRB data, only the TRB chain is used.
#'
#' `cdr3_col` must name a column containing TRB CDR3 amino-acid sequences.
#' Missing and empty CDR3 values are skipped. The other input columns are
#' optional: set `v_gene_col`, `patient_col`, `hla_col`, or `counts_col` only
#' when the analysis should use those data. Every selected column must exist in
#' the `ImmunData` dataset.
#'
#' @section Output from `tool$run()`:
#' `tool$run()` returns an [ImmGliphRun]. It records where the results were
#' saved. To read the result tables, first run `results <- run$results()`.
#'
#' [ImmGliphResults] provides:
#'
#' - **`results$groups`:** One row for each receptor assigned to a group. The
#'   `immgliph_group_id` column gives the group. Receptors without a group are
#'   absent. Use `dplyr::collect()` to read the table into R.
#' - **`results$connections`:** Sequence pairs connected by ImmGLIPH and the
#'   reason for each connection.
#' - **`results$cluster_properties`:** One row per group with the properties and
#'   scores calculated by ImmGLIPH.
#' - **`results$mapping`:** The receptor columns that were passed to ImmGLIPH.
#' - **`results$raw_output`:** Path to the original ImmGLIPH `.rds` file.
#' - **`results$annotate(idata)`:** A new `ImmunData` dataset with an
#'   `immgliph_group_id` receptor column. Receptors without a group have `NA`.
#'
#' @param method ImmGLIPH analysis preset: `"gliph2"`, `"gliph1"`, or
#'   `"custom"`.
#' @param reference Name of the reference repertoire included with ImmGLIPH.
#' @param cdr3_col Name of the column containing TRB CDR3 amino-acid
#'   sequences.
#' @param v_gene_col Optional name of the TRB V-gene column.
#' @param patient_col Optional name of the donor or sample column.
#' @param hla_col Optional name of the HLA column.
#' @param counts_col Optional name of the clone-count column.
#' @param n_cores Positive integer number of CPU cores used by ImmGLIPH.
#' @param seed Optional integer random seed for a reproducible ImmGLIPH run.
#' @param run_args Additional named arguments passed to `immGLIPH::runGLIPH()`.
#'   Use this only for advanced ImmGLIPH settings. Arguments controlled by
#'   immundapter cannot be changed here.
#' @param config An [ImmGliphToolConfig] created by `conf_immgliph()`.
#'
#' @return [conf_immgliph()] returns a read-only [ImmGliphToolConfig].
#'   [tool_immgliph()] returns a reusable [ImmGliphTool].
#'
#' @examplesIf identical(Sys.getenv("IMD_RUN_EXAMPLES"), "true")
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
#'   output_folder = tempfile("immgliph-immundata-")
#' )
#' config <- conf_immgliph(
#'   v_gene_col = "v_call",
#'   n_cores = 4L,
#'   seed = 1L
#' )
#' tool <- tool_immgliph(config)
#' tool$help()
#'
#' if (tool$check()) {
#'   run <- tool$run(idata, output_dir = tempfile("immgliph-run-"))
#'   results <- run$results()
#'
#'   groups <- dplyr::collect(results$groups)
#'   group_scores <- dplyr::collect(results$cluster_properties)
#'   grouped_idata <- results$annotate(idata)
#' }
#' @name immgliph
#' @seealso
#'   [ImmGLIPH manual](https://bioc.r-universe.dev/immGLIPH/doc/manual.html#runGLIPH)
#' @export
conf_immgliph <- function(method = "gliph2",
                          reference = "human_v2.0_CD48",
                          cdr3_col = "junction_aa",
                          v_gene_col = NULL,
                          patient_col = NULL,
                          hla_col = NULL,
                          counts_col = NULL,
                          n_cores = 1L,
                          seed = NULL,
                          run_args = list()) {
  ImmGliphToolConfig$new(list(
    method = method,
    reference = reference,
    cdr3_col = cdr3_col,
    v_gene_col = v_gene_col,
    patient_col = patient_col,
    hla_col = hla_col,
    counts_col = counts_col,
    n_cores = n_cores,
    seed = seed,
    run_args = run_args
  ))
}


#' @rdname immgliph
#' @export
tool_immgliph <- function(config = conf_immgliph()) {
  ImmGliphTool$new(config)
}


#' @rdname immgliph
#'
#' @description
#' `ImmGliphToolConfig` stores the validated settings for an ImmGLIPH analysis.
#' It contains no repertoire data and no results. Create it with
#' [conf_immgliph()] and pass it to [tool_immgliph()].
#'
#' @field method ImmGLIPH preset: `"gliph2"`, `"gliph1"`, or `"custom"`.
#' @field reference Built-in ImmGLIPH reference repertoire name.
#' @field cdr3_col Column containing CDR3 amino-acid sequences.
#' @field v_gene_col Optional V-gene column.
#' @field patient_col Optional donor or sample column.
#' @field hla_col Optional HLA column.
#' @field counts_col Optional clone-count column.
#' @field n_cores Number of cores passed to ImmGLIPH.
#' @field seed Optional random seed used around the ImmGLIPH call.
#' @field run_args Additional named arguments for `immGLIPH::runGLIPH()`.
ImmGliphToolConfig <- R6::R6Class(
  classname = "ImmGliphToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    method = NULL,
    reference = NULL,
    cdr3_col = NULL,
    v_gene_col = NULL,
    patient_col = NULL,
    hla_col = NULL,
    counts_col = NULL,
    n_cores = NULL,
    seed = NULL,
    run_args = NULL,
    #' @description
    #' Create checked ImmGLIPH settings. Users should call
    #' [conf_immgliph()] instead.
    #' @param values A named list of ImmGLIPH settings.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_choice(
        values$method,
        c("gliph2", "gliph1", "custom")
      )
      checkmate::assert_string(values$reference)
      checkmate::assert_string(values$cdr3_col)
      checkmate::assert_string(values$v_gene_col, null.ok = TRUE)
      checkmate::assert_string(values$patient_col, null.ok = TRUE)
      checkmate::assert_string(values$hla_col, null.ok = TRUE)
      checkmate::assert_string(values$counts_col, null.ok = TRUE)
      checkmate::assert_int(values$n_cores, lower = 1L)
      checkmate::assert_int(values$seed, null.ok = TRUE)
      checkmate::assert_list(values$run_args, names = "unique")

      super$initialize(values)
    }
  )
)


#' Run ImmGLIPH analyses
#'
#' `ImmGliphTool` uses one [ImmGliphToolConfig] as its default settings. Create
#' it with [tool_immgliph()]. See [immgliph] for the required `ImmunData` input
#' and available output.
#'
#' @section Methods:
#'
#' - **`$help()`:** Print ImmGLIPH installation instructions.
#' - **`$check()`:** Report whether the immGLIPH R package is available.
#' - **`$run(input, config = NULL, output_dir = NULL)`:** Run ImmGLIPH on an
#'   `ImmunData` dataset and return an [ImmGliphRun].
#' - **`$citation()`:** Reserved for future citation support.
#'
#' @keywords internal
ImmGliphTool <- R6::R6Class(
  classname = "ImmGliphTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an ImmGLIPH tool. Users should call [tool_immgliph()] instead.
    #' @param config An [ImmGliphToolConfig] returned by [conf_immgliph()].
    initialize = function(config = conf_immgliph()) {
      checkmate::assert_r6(config, classes = "ImmGliphToolConfig")
      super$initialize("immgliph", config, ImmGliphRun)
    }
  ),
  private = list(
    help_impl = function() {
      command <- paste0(
        'install.packages("immGLIPH", repos = c(',
        '"https://bioc.r-universe.dev", ',
        '"https://cloud.r-project.org"))'
      )
      documentation <- paste0(
        "https://bioc.r-universe.dev/immGLIPH/",
        "doc/manual.html#runGLIPH"
      )

      cli::cli_inform(c("i" = "Install {.pkg immGLIPH}:"))
      div <- cli::cli_div(theme = list(
        ol = list("margin-left" = 2),
        ul = list("margin-left" = 2)
      ))
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Install the R package on macOS, Linux, or Windows.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {command}}")
      cli::cli_end(commands)
      cli::cli_end(steps)
      cli::cli_end(div)

      cli::cli_inform(c("i" = "Documentation:"))
      div <- cli::cli_div(theme = list(
        ul = list("margin-left" = 2)
      ))
      links <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.url {documentation}}")
      cli::cli_end(links)
      cli::cli_end(div)
    },
    check_impl = function() {
      runtime <- immgliph_runtime_info()

      if (runtime$installed) {
        cli::cli_inform(c(
          "v" = "{.pkg immGLIPH} {runtime$version} is installed.",
          ">" = "Location: {.path {runtime$path}}"
        ))
      } else {
        cli::cli_inform(c(
          "x" = "{.pkg immGLIPH} is not installed.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      runtime$installed
    },
    citation_impl = function() {
      rlang::abort("ImmGLIPH `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      if (!inherits(input, "ImmunData")) {
        rlang::abort("Pass an `ImmunData` dataset to ImmGLIPH.")
      }
      checkmate::assert_r6(input, classes = "ImmunData")
      input
    },
    run_impl = function(input, config, run_dir) {
      if (!immgliph_runtime_info()$installed) {
        private$help_impl()
        rlang::abort("Install immGLIPH before running this tool.")
      }

      prepared <- immgliph_prepare_immundata(input, config)
      mapping <- prepared$mapping

      input <- mapping |>
        dplyr::select(dplyr::all_of(prepared$input_columns)) |>
        dplyr::distinct() |>
        dplyr::collect()

      if (nrow(input) == 0L) {
        rlang::abort("No TRB CDR3 amino-acid sequences are available.")
      }

      input_path <- fs::path(run_dir, "immgliph-input.tsv")
      mapping_path <- fs::path(run_dir, "immgliph-mapping.parquet")
      result_path <- fs::path(run_dir, "immgliph-result.rds")
      groups_path <- fs::path(run_dir, "immgliph-groups.parquet")
      connections_path <- fs::path(run_dir, "immgliph-connections.parquet")
      properties_path <- fs::path(
        run_dir,
        "immgliph-cluster-properties.parquet"
      )

      utils::write.table(
        input,
        file = input_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )
      duckplyr::compute_parquet(mapping, mapping_path)

      arguments <- utils::modifyList(
        config$run_args,
        list(
          cdr3_sequences = input,
          method = config$method,
          chains = "TRB",
          result_folder = run_dir,
          refdb_beta = config$reference,
          n_cores = config$n_cores,
          verbose = FALSE
        )
      )

      result <- if (is.null(config$seed)) {
        immgliph_call(arguments)
      } else {
        withr::with_seed(
          config$seed,
          immgliph_call(arguments)
        )
      }
      saveRDS(result, result_path)

      groups <- immgliph_build_groups(
        result = result,
        mapping = mapping,
        cdr3_col = config$cdr3_col
      )
      duckplyr::compute_parquet(groups, groups_path)

      connections <- result$connections
      if (!is.data.frame(connections)) {
        rlang::abort(
          "ImmGLIPH output `connections` must be a data frame."
        )
      }
      connections <- duckplyr::as_duckdb_tibble(connections)
      duckplyr::compute_parquet(connections, connections_path)

      properties <- result$cluster_properties
      if (!is.data.frame(properties)) {
        rlang::abort(
          "ImmGLIPH output `cluster_properties` must be a data frame."
        )
      }
      properties <- duckplyr::as_duckdb_tibble(properties)
      duckplyr::compute_parquet(properties, properties_path)

      list(
        result = result_path,
        input = input_path,
        mapping = mapping_path,
        groups = groups_path,
        connections = connections_path,
        cluster_properties = properties_path
      )
    }
  )
)
