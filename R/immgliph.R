#' Configure or create an ImmGLIPH adapter
#'
#' `conf_immgliph()` validates ImmGLIPH parameters. `tool_immgliph()` creates
#' the adapter from that configuration.
#'
#' @section Development status:
#' ImmGLIPH configuration, execution, result parsing, and annotation are
#' available. Citation metadata are planned for the next version.
#'
#' @param method ImmGLIPH preset: `"gliph2"`, `"gliph1"`, or `"custom"`.
#' @param reference Built-in ImmGLIPH reference repertoire name.
#' @param cdr3_col Column containing CDR3 amino-acid sequences.
#' @param v_gene_col Optional V-gene column.
#' @param patient_col Optional donor or sample column.
#' @param hla_col Optional HLA column.
#' @param counts_col Optional clone-count column.
#' @param n_cores Number of cores passed to ImmGLIPH.
#' @param seed Optional random seed used around the ImmGLIPH call.
#' @param run_args Additional named arguments passed to `immGLIPH::runGLIPH()`.
#' @param config An [ImmGliphToolConfig] created by `conf_immgliph()`.
#'
#' @return `conf_immgliph()` returns an ImmGLIPH configuration.
#'   `tool_immgliph()` returns an ImmGLIPH adapter.
#'
#' @seealso \url{https://bioc.r-universe.dev/immGLIPH/doc/manual.html#runGLIPH}
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


#' @rdname conf_immgliph
#' @export
tool_immgliph <- function(config = conf_immgliph()) {
  ImmGliphTool$new(config)
}


#' ImmGLIPH tool configuration
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
#'
#' @keywords internal
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
    #' Create a validated ImmGLIPH configuration. Users should call
    #' [conf_immgliph()] instead.
    #' @param values A named list of ImmGLIPH configuration values.
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


#' ImmGLIPH adapter
#'
#' @keywords internal
ImmGliphTool <- R6::R6Class(
  classname = "ImmGliphTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an ImmGLIPH adapter. Users should call [tool_immgliph()] instead.
    #' @param config An [ImmGliphToolConfig] returned by [conf_immgliph()].
    initialize = function(config = conf_immgliph()) {
      checkmate::assert_r6(config, classes = "ImmGliphToolConfig")
      super$initialize("immgliph", config)
    }
  ),
  private = list(
    help_impl = function() {
      command <- 'pak::pak("BorchLab/immGLIPH")'
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
      installed <- requireNamespace("immGLIPH", quietly = TRUE) &&
        exists("runGLIPH", envir = asNamespace("immGLIPH"), inherits = FALSE)

      if (installed) {
        version <- as.character(utils::packageVersion("immGLIPH"))
        path <- system.file(package = "immGLIPH")
        cli::cli_inform(c(
          "v" = "{.pkg immGLIPH} {version} is installed.",
          ">" = "Location: {.path {path}}"
        ))
      } else {
        cli::cli_inform(c(
          "x" = "{.pkg immGLIPH} is not installed.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      installed
    },
    citation_impl = function() {
      rlang::abort("ImmGLIPH `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      checkmate::assert_r6(input, classes = "ImmunData")
      input
    },
    run_impl = function(input, config, run_dir) {
      if (!requireNamespace("immGLIPH", quietly = TRUE)) {
        private$help_impl()
        rlang::abort("immGLIPH is required to run this adapter.")
      }

      idata <- input
      receptor_col <- immundata::imd_schema("receptor")
      locus_col <- immundata::imd_schema("locus")
      input_columns <- c(
        CDR3b = config$cdr3_col,
        TRBV = config$v_gene_col,
        patient = config$patient_col,
        HLA = config$hla_col,
        counts = config$counts_col
      )
      selected_columns <- unique(c(receptor_col, unname(input_columns)))

      mapping <- idata$annotations |>
        dplyr::filter(!!rlang::sym(locus_col) == "TRB") |>
        dplyr::select(dplyr::all_of(selected_columns)) |>
        dplyr::filter(!is.na(!!rlang::sym(config$cdr3_col))) |>
        dplyr::filter(!!rlang::sym(config$cdr3_col) != "") |>
        dplyr::distinct()

      input <- mapping |>
        dplyr::select(dplyr::all_of(input_columns)) |>
        dplyr::distinct() |>
        dplyr::collect()

      if (nrow(input) == 0L) {
        rlang::abort("No TRB CDR3 amino-acid sequences are available.")
      }

      input_path <- fs::path(run_dir, "immgliph-input.tsv")
      mapping_path <- fs::path(run_dir, "immgliph-mapping.parquet")
      result_path <- fs::path(run_dir, "immgliph-result.rds")

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
        do.call(getExportedValue("immGLIPH", "runGLIPH"), arguments)
      } else {
        withr::with_seed(
          config$seed,
          do.call(getExportedValue("immGLIPH", "runGLIPH"), arguments)
        )
      }
      saveRDS(result, result_path)

      list(
        result = result_path,
        input = input_path,
        mapping = mapping_path
      )
    },
    results_impl = function(run) {
      checkmate::assert_file_exists(run$artifacts$result)
      checkmate::assert_file_exists(run$artifacts$mapping)

      raw <- readRDS(run$artifacts$result)
      clusters <- raw$cluster_list

      if (length(clusters) == 0L) {
        membership <- data.frame(
          CDR3b = character(),
          immgliph_group_id = character()
        )
      } else {
        membership <- do.call(
          rbind,
          Map(
            function(cluster, group_id) {
              data.frame(
                CDR3b = as.character(cluster$CDR3b),
                immgliph_group_id = group_id
              )
            },
            clusters,
            names(clusters)
          )
        )
      }

      membership <- duckplyr::as_duckdb_tibble(membership)
      mapping <- duckplyr::read_parquet_duckdb(run$artifacts$mapping)
      receptor_col <- immundata::imd_schema("receptor")

      mapping |>
        dplyr::select(dplyr::all_of(c(receptor_col, run$config$cdr3_col))) |>
        dplyr::rename(CDR3b = dplyr::all_of(run$config$cdr3_col)) |>
        dplyr::inner_join(membership, by = "CDR3b") |>
        dplyr::select(
          dplyr::all_of(c(receptor_col, "immgliph_group_id"))
        ) |>
        dplyr::distinct()
    },
    annotate_impl = function(idata, results, ...) {
      immundata::annotate_receptors(
        idata = idata,
        annotations = results,
        annot_col = immundata::imd_schema("receptor"),
        ...
      )
    }
  )
)
