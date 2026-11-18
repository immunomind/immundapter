#' @title Compare and group T-cell receptors with tcrdist3
#'
#' @description
#' tcrdist3 measures how similar T-cell receptors are from their CDR3
#' amino-acid sequences and V/J genes. A smaller distance means that two
#' receptors are more similar. The tool can return the distances alone or use
#' them to place similar receptors into groups.
#'
#' The usual workflow is:
#'
#' ```r
#' tool <- tool_tcrdist3()
#' run <- tool$run(idata)
#' results <- run$results()
#' clusters <- dplyr::collect(results$clusters)
#' ```
#'
#' @details
#' tcrdist3 compares receptor sequences; it does not prove that receptors bind
#' the same antigen. Treat the distance and cluster threshold as analysis
#' choices that should be checked for the biological question and dataset.
#'
#' Every selected chain needs a CDR3 amino-acid sequence, a V-gene call, and a
#' J-gene call. The defaults expect human TRB data stored in `junction_aa`,
#' `v_call`, and `j_call`.
#'
#' @section Choose distances or clusters:
#'
#' Use `mode = "cluster"`, the default, to calculate distances and group
#' receptors. Two receptors are connected when their distance is at most
#' `radius`. All receptors connected directly or through other receptors become
#' one cluster. Consequently, two receptors in the same cluster can have a
#' direct distance greater than `radius`.
#'
#' Use `mode = "distance"` when only the pairwise distance matrix is needed.
#' Distance-only results do not contain clusters and cannot be added to
#' `ImmunData` with `$annotate()`.
#'
#' @section Choose the receptor chains:
#'
#' `chains` selects the chains used in the distance:
#'
#' * `chains = "beta"`, the default, compares TRB only. It can select the beta
#'   chain from paired TRA/TRB `ImmunData`.
#' * `chains = c("alpha", "beta")` compares both sides of a paired receptor.
#'   tcrdist3 calculates one distance per chain and adds them.
#'
#' Alpha, beta, gamma, and delta chains are supported. `organism` must match the
#' gene names in the data. Because paired-chain distances are sums, a `radius`
#' chosen for one chain should not be assumed to work for two chains.
#'
#' For paired `ImmunData`, the receptor schema tells the tool which row belongs
#' to each chain, so both chains can use the same feature names:
#'
#' ```r
#' conf_tcrdist3(
#'   chains = c("alpha", "beta"),
#'   cdr3_cols = c(alpha = "junction_aa", beta = "junction_aa"),
#'   v_gene_cols = c(alpha = "v_call", beta = "v_call"),
#'   j_gene_cols = c(alpha = "j_call", beta = "j_call")
#' )
#' ```
#'
#' @section Choose the input:
#'
#' `tool$run()` accepts:
#'
#' * **An `ImmunData` dataset.** The selected columns must be receptor features.
#'   The tool prepares one row per biological receptor. Cluster results can be
#'   added back with `results$annotate(idata)`.
#' * **A CSV, TSV, or TXT file.** Each row must represent one receptor. Paired
#'   receptors need separate CDR3, V-gene, and J-gene columns for each chain.
#'
#' Rows missing a required value for a selected chain are excluded. For text
#' input, `imd_receptor_id` is used as the identifier when present, followed by
#' `cell_id`; otherwise, row numbers are used.
#'
#' @section Run the tool:
#'
#' `tool_tcrdist3()` creates a reusable tool. Its main methods are:
#'
#' * `tool$help()`: print commands for installing the Python environment;
#' * `tool$check()`: check that Python and tcrdist3 are ready;
#' * `tool$run(input, config = NULL, output_dir = NULL)`: run one analysis.
#'
#' The tool keeps its default settings in `tool$config`. Supplying `config` to
#' `$run()` replaces those settings for that run only. The completed run stores
#' its output directory in `$path` and its saved files in `$artifacts`.
#'
#' @section Read the results:
#'
#' Call `run$results()` after a successful run. Both modes provide:
#'
#' * `results$distances`: the path to the complete NumPy distance matrix;
#' * `results$mapping`: the link between input receptors and rows of that
#'   matrix.
#'
#' Cluster mode also provides:
#'
#' * `results$clusters`: one cluster number and cluster size per receptor;
#' * `results$edges`: receptor pairs whose distance is at most `radius`;
#' * `results$annotate(idata)`: a new `ImmunData` dataset containing
#'   `tcrdist_cluster_id` and `tcrdist_cluster_size`.
#'
#' Result tables stay on disk and are evaluated with duckplyr. Use
#' `dplyr::collect()` when a table is small enough to bring into R memory. The
#' distance matrix is a NumPy `.npy` file, not an R table.
#'
#' @section Calculation details:
#'
#' With `deduplicate = TRUE`, identical selected receptors are calculated once
#' and then matched back to every original receptor. `results$mapping` records
#' this relationship. Set `store_all_cdr = FALSE` when separate
#' region-specific distance matrices are not needed.
#'
#' @section Python setup:
#'
#' By default, immundapter uses a dedicated Python environment at
#' `environments/tcrdist3` under `tools::R_user_dir("immundapter", "data")`.
#' Call `tool$help()` to print the installation commands and `tool$check()`
#' to verify the environment.
#' `$help()` prints instructions; it does not install software.
#'
#' Set `IMMUNDAPTER_DATA` before creating the configuration to use another data
#' directory, or supply `python` to use a specific Python executable.
#'
#' @param mode `"cluster"` groups similar receptors and returns their distances.
#'   `"distance"` returns only the distance matrix and its input mapping.
#' @param organism Species used to interpret V/J genes: `"human"` or `"mouse"`.
#' @param chains One or two chains selected from `"alpha"`, `"beta"`, `"gamma"`,
#'   and `"delta"`. With two chains, their distances are added.
#' @param cdr3_cols Named vector giving the CDR3 amino-acid column for each
#'   selected chain, for example `c(beta = "junction_aa")`.
#' @param v_gene_cols Named vector giving the V-gene column for each selected
#'   chain. Gene names must be recognized by the tcrdist3 database.
#' @param j_gene_cols Named vector giving the J-gene column for each selected
#'   chain.
#' @param db_file Germline reference supplied to tcrdist3. The default reference
#'   is included with tcrdist3.
#' @param store_all_cdr Whether to save separate distance matrices for
#'   individual receptor regions in addition to the combined matrices.
#' @param deduplicate Whether to calculate identical receptors once. Results
#'   are still reported for every original receptor.
#' @param cpus Number of CPU cores used by tcrdist3.
#' @param radius Largest distance that directly connects two receptors during
#'   clustering. It is used only when `mode = "cluster"`.
#' @param python Python executable name or path. `NULL` selects the managed
#'   tcrdist3 environment described above.
#' @param run_args Additional named arguments passed to `TCRrep`. This is for
#'   advanced settings. Values must be valid JSON values. Arguments controlled
#'   by immundapter cannot be changed here.
#'
#' @param config Settings created by [conf_tcrdist3()].
#'
#' @return [conf_tcrdist3()] returns read-only tcrdist3 settings.
#'   [tool_tcrdist3()] returns a reusable tcrdist3 tool.
#'
#' @examplesIf identical(Sys.getenv("IMD_RUN_EXAMPLES"), "true")
#' # Group similar TRB receptors in the bundled example.
#' config <- conf_tcrdist3(
#'   store_all_cdr = FALSE,
#'   cpus = 1L,
#'   radius = 24L
#' )
#' tool <- tool_tcrdist3(config)
#' stopifnot(tool$check())
#'
#' input_file <- system.file(
#'   "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
#'   package = "immundapter"
#' )
#' run <- tool$run(
#'   input_file,
#'   output_dir = tempfile("tcrdist3-file-")
#' )
#' results <- run$results()
#' clusters <- dplyr::collect(results$clusters)
#'
#' # Add the cluster assignments to ImmunData.
#' snapshot <- system.file(
#'   "extdata", "tcrdist3", "immundata-tcr-cluster-simple",
#'   package = "immundapter"
#' )
#' idata <- immundata::read_immundata(snapshot, verbose = FALSE)
#' idata_run <- tool$run(
#'   idata,
#'   output_dir = tempfile("tcrdist3-immundata-")
#' )
#' annotated_idata <- idata_run$results()$annotate(idata)
#'
#' # Calculate distances without creating clusters.
#' distance_tool <- tool_tcrdist3(conf_tcrdist3(
#'   mode = "distance", store_all_cdr = FALSE
#' ))
#' distance_run <- distance_tool$run(
#'   input_file,
#'   output_dir = tempfile("tcrdist3-distance-")
#' )
#' distance_file <- distance_run$results()$distances
#' @name tool_tcrdist3
#' @seealso [immundata::ImmunData],
#'   [tcrdist3 input](https://tcrdist3.readthedocs.io/en/latest/inputs.html),
#'   [tcrdist3 documentation](https://tcrdist3.readthedocs.io/en/latest/)
#' @export
conf_tcrdist3 <- function(mode = c("cluster", "distance"),
                          organism = "human",
                          chains = "beta",
                          cdr3_cols = c(beta = "junction_aa"),
                          v_gene_cols = c(beta = "v_call"),
                          j_gene_cols = c(beta = "j_call"),
                          db_file = "alphabeta_gammadelta_db.tsv",
                          store_all_cdr = TRUE,
                          deduplicate = TRUE,
                          cpus = 1L,
                          radius = 24L,
                          python = NULL,
                          run_args = list()) {
  mode <- match.arg(mode)
  if (is.null(python)) {
    python <- tcrdist3_managed_python()
  }

  Tcrdist3ToolConfig$new(list(
    mode = mode,
    organism = organism,
    chains = chains,
    cdr3_cols = cdr3_cols,
    v_gene_cols = v_gene_cols,
    j_gene_cols = j_gene_cols,
    db_file = db_file,
    store_all_cdr = store_all_cdr,
    deduplicate = deduplicate,
    cpus = cpus,
    radius = radius,
    python = python,
    run_args = run_args
  ))
}


#' @rdname tool_tcrdist3
#' @export
tool_tcrdist3 <- function(config = conf_tcrdist3()) {
  Tcrdist3Tool$new(config)
}


#' @noRd
Tcrdist3ToolConfig <- R6::R6Class(
  classname = "Tcrdist3ToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    mode = NULL,
    organism = NULL,
    chains = NULL,
    cdr3_cols = NULL,
    v_gene_cols = NULL,
    j_gene_cols = NULL,
    db_file = NULL,
    store_all_cdr = NULL,
    deduplicate = NULL,
    cpus = NULL,
    radius = NULL,
    python = NULL,
    run_args = NULL,
    #' @description
    #' Create checked tcrdist3 settings. Users should call
    #' [conf_tcrdist3()] instead.
    #' @param values Named list of tcrdist3 settings.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_choice(values$mode, c("cluster", "distance"))
      checkmate::assert_choice(values$organism, c("human", "mouse"))
      checkmate::assert_character(
        values$chains,
        min.len = 1L,
        max.len = 2L,
        unique = TRUE
      )
      checkmate::assert_subset(
        values$chains,
        c("alpha", "beta", "gamma", "delta")
      )
      checkmate::assert_character(values$cdr3_cols, names = "unique")
      checkmate::assert_character(values$v_gene_cols, names = "unique")
      checkmate::assert_character(values$j_gene_cols, names = "unique")
      checkmate::assert_names(
        names(values$cdr3_cols),
        permutation.of = values$chains
      )
      checkmate::assert_names(
        names(values$v_gene_cols),
        permutation.of = values$chains
      )
      checkmate::assert_names(
        names(values$j_gene_cols),
        permutation.of = values$chains
      )
      checkmate::assert_string(values$db_file)
      checkmate::assert_flag(values$store_all_cdr)
      checkmate::assert_flag(values$deduplicate)
      checkmate::assert_int(values$cpus, lower = 1L)
      checkmate::assert_int(values$radius, lower = 0L)
      checkmate::assert_string(values$python)
      checkmate::assert_list(values$run_args, names = "unique")

      super$initialize(values)
    }
  )
)


tcrdist3_managed_python <- function() {
  executable <- if (.Platform$OS.type == "windows") {
    fs::path("Scripts", "python.exe")
  } else {
    fs::path("bin", "python")
  }
  immundapter_data_path("environments", "tcrdist3", executable, which = "data")
}

#' @noRd
Tcrdist3Tool <- R6::R6Class(
  classname = "Tcrdist3Tool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a tcrdist3 tool. Users should call [tool_tcrdist3()] instead.
    #' @param config Settings returned by [conf_tcrdist3()].
    initialize = function(config = conf_tcrdist3()) {
      checkmate::assert_r6(config, classes = "Tcrdist3ToolConfig")
      super$initialize("tcrdist3", config, Tcrdist3Run)
    }
  ),
  private = list(
    help_impl = function() {
      python <- Sys.which(private$.config$python)
      uv <- Sys.which("uv")
      configured_python <- private$.config$python

      if (.Platform$OS.type == "windows") {
        install_uv <- paste0(
          'powershell -ExecutionPolicy ByPass -c "',
          'irm https://astral.sh/uv/install.ps1 | iex"'
        )
        shell <- "Windows PowerShell"
        quote_path <- function(path) {
          paste0("'", gsub("'", "''", path, fixed = TRUE), "'")
        }
      } else {
        install_uv <- "curl -LsSf https://astral.sh/uv/install.sh | sh"
        shell <- "a macOS or Linux terminal"
        quote_path <- shQuote
      }
      managed_python <- tcrdist3_managed_python()
      environment_dir <- fs::path_dir(fs::path_dir(managed_python))
      create_env <- paste("uv venv", quote_path(environment_dir))
      install_tcrdist3 <- paste(
        "uv pip install --python", quote_path(managed_python), "tcrdist3==0.3"
      )
      uv_documentation <- paste0(
        "https://docs.astral.sh/uv/getting-started/",
        "installation/"
      )
      tcrdist3_documentation <- paste0(
        "https://tcrdist3.readthedocs.io/en/latest/"
      )

      cli::cli_inform(c(
        "i" = paste0(
          "Create an isolated Python environment for ",
          "{.pkg tcrdist3}:"
        )
      ))
      cli::cli_inform(c(
        if (nzchar(uv)) {
          uv_version <- system2(
            uv,
            "--version",
            stdout = TRUE,
            stderr = TRUE
          )[[1L]]
          c(
            "v" = "{uv_version}",
            ">" = "uv executable: {.path {uv}}"
          )
        } else {
          c("!" = "{.pkg uv} is not currently installed.")
        },
        if (nzchar(python)) {
          python_version <- system2(
            python,
            "--version",
            stdout = TRUE,
            stderr = TRUE
          )[[1L]]
          c(
            "v" = "{python_version}",
            ">" = "Python executable: {.path {python}}"
          )
        } else {
          c(
            "!" = paste0(
              "The configured Python does not exist yet: ",
              "{.path {configured_python}}"
            )
          )
        }
      ))

      div <- cli::cli_div(theme = list(
        ol = list("margin-left" = 2),
        ul = list("margin-left" = 2)
      ))
      cli::cli_inform("Run these commands in {shell}.")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Install {.pkg uv}.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {install_uv}}")
      cli::cli_end(commands)
      cli::cli_li("Create the Python environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {create_env}}")
      cli::cli_end(commands)
      cli::cli_li("Install {.pkg tcrdist3} into the environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {install_tcrdist3}}")
      cli::cli_end(commands)
      cli::cli_end(steps)
      cli::cli_end(div)

      cli::cli_inform(c(
        "i" = paste0(
          "Set {.envvar IMMUNDAPTER_DATA} before creating the tool ",
          "to use another data directory."
        )
      ))
      cli::cli_inform(c("i" = "Documentation:"))
      div <- cli::cli_div(theme = list(
        ul = list("margin-left" = 2)
      ))
      links <- cli::cli_ul(.close = FALSE)
      cli::cli_li("uv: {.url {uv_documentation}}")
      cli::cli_li("tcrdist3: {.url {tcrdist3_documentation}}")
      cli::cli_end(links)
      cli::cli_end(div)
    },
    check_impl = function() {
      python <- Sys.which(private$.config$python)
      python_version <- NULL
      installed <- FALSE

      if (nzchar(python)) {
        python_version <- system2(
          python,
          "--version",
          stdout = TRUE,
          stderr = TRUE
        )[[1L]]
        status <- suppressWarnings(system2(
          python,
          c("-c", shQuote("import tcrdist")),
          stdout = FALSE,
          stderr = FALSE
        ))
        installed <- identical(status, 0L)
      }

      if (installed) {
        cli::cli_inform(c(
          "v" = "{.pkg tcrdist3} is installed.",
          "v" = "{python_version}",
          ">" = "Python executable: {.path {python}}"
        ))
      } else if (nzchar(python)) {
        cli::cli_inform(c(
          "v" = "{python_version}",
          ">" = "Python executable: {.path {python}}",
          "x" = "{.pkg tcrdist3} is not installed in that environment.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      } else {
        executable <- private$.config$python
        cli::cli_inform(c(
          "x" = "Python executable {.code {executable}} was not found.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      installed
    },
    citation_impl = function() {
      rlang::abort("tcrdist3 `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      if (inherits(input, "ImmunData")) {
        checkmate::assert_r6(input, classes = "ImmunData")
        return(input)
      }
      checkmate::assert_string(input, min.chars = 1L)
      checkmate::assert_file_exists(input)
      fs::path_abs(input)
    },
    run_impl = function(input, config, run_dir) {
      python <- Sys.which(config$python)
      if (!nzchar(python)) {
        rlang::abort(sprintf(
          "The configured tcrdist3 Python executable was not found: %s.",
          config$python
        ))
      }

      import_status <- suppressWarnings(system2(
        python,
        c("-c", shQuote("import tcrdist")),
        stdout = FALSE,
        stderr = FALSE
      ))
      if (!identical(import_status, 0L)) {
        rlang::abort(
          "tcrdist3 is not installed in the configured Python environment."
        )
      }

      prepared <- tcrdist3_prepare_input(input, config)
      prepared_df <- as.data.frame(dplyr::collect(prepared))
      if (nrow(prepared_df) == 0L) {
        rlang::abort("No complete TCRs are available for tcrdist3.")
      }
      prepared_df$input_id <- as.character(prepared_df$input_id)

      input_path <- fs::path(run_dir, "tcrdist3-input.tsv")
      config_path <- fs::path(run_dir, "tcrdist3-config.json")
      clusters_path <- fs::path(run_dir, "tcrdist3-clusters.tsv")
      mapping_path <- fs::path(run_dir, "tcrdist3-mapping.tsv")
      edges_path <- fs::path(run_dir, "tcrdist3-edges.tsv")
      distances_path <- fs::path(run_dir, "tcrdist3-distances.npy")
      stdout_path <- fs::path(run_dir, "stdout.log")
      stderr_path <- fs::path(run_dir, "stderr.log")

      utils::write.table(
        prepared_df,
        file = input_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )
      jsonlite::write_json(
        list(
          mode = config$mode,
          organism = config$organism,
          chains = config$chains,
          db_file = config$db_file,
          store_all_cdr = config$store_all_cdr,
          deduplicate = config$deduplicate,
          cpus = config$cpus,
          radius = config$radius,
          run_args = config$run_args
        ),
        path = config_path,
        auto_unbox = TRUE,
        null = "null",
        pretty = TRUE
      )

      runner <- system.file(
        "python",
        "tcrdist3_runner.py",
        package = "immundapter"
      )
      if (!nzchar(runner)) {
        rlang::abort("The bundled tcrdist3 runner was not found.")
      }

      status <- system2(
        python,
        c(
          shQuote(runner),
          "--input", shQuote(input_path),
          "--config", shQuote(config_path),
          "--clusters", shQuote(clusters_path),
          "--mapping", shQuote(mapping_path),
          "--edges", shQuote(edges_path),
          "--distances", shQuote(distances_path)
        ),
        stdout = stdout_path,
        stderr = stderr_path
      )
      if (!identical(status, 0L)) {
        details <- if (fs::file_exists(stderr_path)) {
          paste(tail(readLines(stderr_path, warn = FALSE), 20L), collapse = "\n")
        } else {
          "No stderr log was produced."
        }
        rlang::abort(paste(
          "tcrdist3 execution failed.",
          details,
          sep = "\n"
        ))
      }

      artifacts <- list(
        input = input_path,
        config = config_path,
        mapping = mapping_path,
        distances = distances_path,
        stdout = stdout_path,
        stderr = stderr_path
      )
      if (identical(config$mode, "cluster")) {
        artifacts$clusters <- clusters_path
        artifacts$edges <- edges_path
      }
      artifacts
    }
  )
)
