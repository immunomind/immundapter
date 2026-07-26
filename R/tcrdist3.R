#' Cluster T-cell receptors with tcrdist3
#'
#' tcrdist3 measures similarity between T-cell receptors from their CDR3
#' amino-acid sequences and V/J genes. This adapter prepares R data for the
#' Python `tcrdist3` package, calculates receptor distances, turns the distance
#' matrix into clusters, and brings the cluster assignments back into R.
#'
#' @section The objects and how they connect:
#' A typical analysis uses four objects:
#' \describe{
#'   \item{Configuration}{[conf_tcrdist3()] validates analysis settings and
#'   returns a read-only [Tcrdist3ToolConfig]. The configuration says which
#'   chains contribute to distance, where their CDR3 and V/J columns are, and
#'   how tcrdist3 should run. It contains no repertoire data.}
#'   \item{Tool}{[tool_tcrdist3()] wraps a configuration in a
#'   [Tcrdist3Tool]. The tool can be reused for multiple inputs. It provides
#'   `$help()`, `$check()`, and `$run()`.}
#'   \item{Run}{`tool$run(input)` executes one analysis and returns a
#'   `ToolRun`. The run records the configuration, input summary, run
#'   directory, and paths to every generated artifact. Call `run$results()` to
#'   parse that specific execution.}
#'   \item{Results}{`run$results()` returns [Tcrdist3Results], with lazy
#'   duckplyr `$clusters`, `$mapping`, and `$edges` tables and a `$distances`
#'   artifact path. `results$annotate(idata)` joins clusters back by
#'   `imd_receptor_id` and returns a new annotated `ImmunData`.}
#' }
#'
#' In short:
#'
#' ```
#' config -> tool -> run -> results -> annotated ImmunData
#' ```
#'
#' @section Choosing single- or paired-chain distance:
#' `chains` is passed to `tcrdist.repertoire.TCRrep`. It controls the distance
#' calculation, not the shape of the stored data:
#' \describe{
#'   \item{`chains = "beta"`}{Calculate distance from the beta chain only. This
#'   can also be used with paired TRA/TRB `ImmunData`; the beta side is selected
#'   from each paired receptor.}
#'   \item{`chains = c("alpha", "beta")`}{Calculate alpha and beta distance
#'   matrices and add them before clustering. The input must provide both
#'   chains.}
#' }
#'
#' `ImmunData$schema_receptor` describes the chains available in the data.
#' `config$chains` describes the chains used by this particular analysis. When
#' the schema declares its chains, the adapter checks that the requested chains
#' are compatible with them.
#'
#' @section Input:
#' The recommended input is an [immundata::ImmunData] object. The adapter reads
#' `input$receptors`, which has one row per distinct `imd_receptor_id`.
#' Single-chain receptors use unsuffixed schema feature columns. Paired
#' receptors use the schema loci to associate the `.x` and `.y` columns with
#' alpha, beta, gamma, or delta. The source column names in `cdr3_cols`,
#' `v_gene_cols`, and `j_gene_cols` refer to the unsuffixed schema feature
#' names. For example, paired TRA/TRB `ImmunData` commonly uses:
#'
#' ```
#' cdr3_cols  = c(alpha = "junction_aa", beta = "junction_aa")
#' v_gene_cols = c(alpha = "v_call", beta = "v_call")
#' j_gene_cols = c(alpha = "j_call", beta = "j_call")
#' ```
#'
#' Text files are also supported. CSV files are comma-separated; `.tsv` and
#' `.txt` files are tab-separated. Each row is one original TCR. Paired-chain
#' files must be wide, with separate source columns for every chain. The
#' adapter uses `imd_receptor_id`, then `cell_id`, as `input_id`; if neither is
#' present, it assigns row numbers.
#'
#' @section What happens during a run:
#' The adapter performs these steps:
#' \enumerate{
#'   \item Normalize the selected chains to the column names required by
#'   tcrdist3, such as `cdr3_b_aa`, `v_b_gene`, and `j_b_gene`.
#'   \item If `deduplicate = TRUE`, calculate each identical normalized receptor
#'   only once while retaining an original-to-deduplicated mapping.
#'   \item Run `tcrdist.repertoire.TCRrep` in the configured Python environment.
#'   The adapter has already resolved deduplication, so the runner calls
#'   `TCRrep` with its own deduplication disabled.
#'   \item Calculate one pairwise distance matrix per selected chain and add
#'   the matrices for a multi-chain analysis.
#'   \item Connect receptor pairs whose combined distance is at most `radius`;
#'   each connected component becomes a cluster.
#'   \item Expand cluster assignments, edges, and distances back to the
#'   original input order.
#' }
#'
#' Clustering is transitive: two receptors in one connected component can have
#' a direct distance greater than `radius` when other receptors connect them.
#'
#' @section Managed Python environment:
#' When `python = NULL`, the tool uses an isolated uv environment under
#' `~/immundapter_data/environments/tcrdist3`. On Windows, it uses the
#' corresponding `Scripts/python.exe` path. Set `IMMUNDAPTER_DATA` before
#' creating the configuration to use another root directory. `tool$help()`
#' prints the uv commands that create the environment and install the supported
#' tcrdist3 version.
#'
#' @param organism Host organism: `"human"` or `"mouse"`. It must match the V/J
#'   gene names in the input.
#' @param chains One or two unique chains selected from `"alpha"`, `"beta"`,
#'   `"gamma"`, and `"delta"`. These are passed to `TCRrep` and determine which
#'   chain distances contribute to clustering.
#' @param cdr3_cols Named character vector mapping every value of `chains` to
#'   its source CDR3 amino-acid column. For long-form paired `ImmunData`, two
#'   chains may map to the same source column because `locus` distinguishes
#'   them.
#' @param v_gene_cols Named character vector mapping every value of `chains` to
#'   its source V-gene column. Alleles must be recognized by the configured
#'   tcrdist3 germline database.
#' @param j_gene_cols Named character vector mapping every value of `chains` to
#'   its source J-gene column.
#' @param db_file Germline database passed to
#'   `tcrdist.repertoire.TCRrep`. tcrdist3 includes the default database.
#' @param compute_distances Whether `TCRrep` computes pairwise distances when
#'   it is created. If `FALSE`, the tool computes them in the next step because
#'   clustering requires the complete matrix.
#' @param store_all_cdr Whether `TCRrep` retains the CDR-specific distance
#'   matrices in addition to its per-chain matrices.
#' @param deduplicate Whether identical normalized receptors are calculated
#'   once. Results are expanded back to every original `input_id`.
#' @param cpus Positive integer number of CPUs supplied to `TCRrep`.
#' @param radius Non-negative integer maximum combined tcrdist for an edge in
#'   the clustering graph.
#' @param python Python executable name or path. `NULL` selects the managed
#'   tcrdist3 environment described above.
#' @param run_args Additional named arguments passed to `TCRrep`. The values
#'   must be valid JSON values. Do not include arguments controlled by the
#'   tool, such as `cell_df`, `chains`, or `index_cols`.
#'
#' @param config A [Tcrdist3ToolConfig] returned by [conf_tcrdist3()].
#'
#' @return [conf_tcrdist3()] returns a read-only
#'   [Tcrdist3ToolConfig]. [tool_tcrdist3()] returns a reusable
#'   [Tcrdist3Tool].
#'
#' @examples
#' \dontrun{
#' # 1. Describe the analysis.
#' config <- conf_tcrdist3(cpus = 4L, radius = 24L)
#'
#' # 2. Create a reusable tool.
#' tool <- tool_tcrdist3(config)
#'
#' # One-time environment setup and verification.
#' tool$help()
#' tool$check()
#'
#' # 3. Run the tool on a file.
#' input_file <- system.file(
#'   "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
#'   package = "immundapter"
#' )
#' run <- tool$run(
#'   input_file,
#'   output_dir = file.path(tempdir(), "tcrdist3-file")
#' )
#'
#' # 4. Read the cluster table.
#' results <- run$results()
#' dplyr::collect(results$clusters)
#'
#' # The same configuration can run on ImmunData.
#' snapshot <- system.file(
#'   "extdata", "tcrdist3", "immundata-tcr-cluster-simple",
#'   package = "immundapter"
#' )
#' idata <- immundata::read_immundata(snapshot)
#' idata_run <- tool$run(
#'   idata,
#'   output_dir = file.path(tempdir(), "tcrdist3-immundata")
#' )
#' idata_results <- idata_run$results()
#'
#' # 5. Join receptor-level clusters back to all annotation rows.
#' annotated <- idata_results$annotate(idata)
#'
#' # Build paired TRA/TRB ImmunData from the bundled long-form example.
#' paired_input <- system.file(
#'   "extdata", "tcrdist3", "tcr-cluster-paired-simple.tsv",
#'   package = "immundapter"
#' )
#' paired_idata <- immundata::read_repertoires(
#'   path = paired_input,
#'   schema = immundata::make_receptor_schema(
#'     features = c("v_call", "j_call", "junction_aa"),
#'     chains = c("TRA", "TRB")
#'   ),
#'   barcode_col = "cell_id",
#'   locus_col = "locus",
#'   umi_col = "duplicate_count",
#'   output_folder = file.path(tempdir(), "paired-immundata"),
#'   repertoire_schema = NULL,
#'   preprocess = NULL,
#'   postprocess = NULL,
#'   rename_columns = NULL
#' )
#'
#' # Calculate distance from both sides of each paired receptor.
#' paired_config <- conf_tcrdist3(
#'   chains = c("alpha", "beta"),
#'   cdr3_cols = c(alpha = "junction_aa", beta = "junction_aa"),
#'   v_gene_cols = c(alpha = "v_call", beta = "v_call"),
#'   j_gene_cols = c(alpha = "j_call", beta = "j_call")
#' )
#' paired_tool <- tool_tcrdist3(paired_config)
#' paired_run <- paired_tool$run(
#'   paired_idata,
#'   output_dir = file.path(tempdir(), "tcrdist3-paired")
#' )
#' }
#'
#' @name tcrdist3
#' @seealso [Tcrdist3Tool], [immundata::ImmunData],
#'   \url{https://tcrdist3.readthedocs.io/en/latest/inputs.html},
#'   \url{https://tcrdist3.readthedocs.io/en/latest/}
#' @export
conf_tcrdist3 <- function(organism = "human",
                          chains = "beta",
                          cdr3_cols = c(beta = "junction_aa"),
                          v_gene_cols = c(beta = "v_call"),
                          j_gene_cols = c(beta = "j_call"),
                          db_file = "alphabeta_gammadelta_db.tsv",
                          compute_distances = FALSE,
                          store_all_cdr = TRUE,
                          deduplicate = TRUE,
                          cpus = 1L,
                          radius = 24L,
                          python = NULL,
                          run_args = list()) {
  if (is.null(python)) {
    python <- immundapter_data_path(
      "environments",
      "tcrdist3",
      if (.Platform$OS.type == "windows") {
        fs::path("Scripts", "python.exe")
      } else {
        fs::path("bin", "python")
      }
    )
  }

  Tcrdist3ToolConfig$new(list(
    organism = organism,
    chains = chains,
    cdr3_cols = cdr3_cols,
    v_gene_cols = v_gene_cols,
    j_gene_cols = j_gene_cols,
    db_file = db_file,
    compute_distances = compute_distances,
    store_all_cdr = store_all_cdr,
    deduplicate = deduplicate,
    cpus = cpus,
    radius = radius,
    python = python,
    run_args = run_args
  ))
}


#' @rdname tcrdist3
#' @export
tool_tcrdist3 <- function(config = conf_tcrdist3()) {
  Tcrdist3Tool$new(config)
}


#' @rdname tcrdist3
#'
#' @description
#' `Tcrdist3ToolConfig` is an immutable analysis recipe. It contains validated
#' tcrdist3 settings but no input data or results. Create it with
#' [conf_tcrdist3()], then pass it to [tool_tcrdist3()]. A `ToolRun` stores a
#' plain-list snapshot of the configuration used for that run and exposes its
#' results through `$results()`.
#'
#' @field organism Validated organism supplied to `TCRrep`.
#' @field chains Validated TCR chains.
#' @field cdr3_cols Named chain-to-CDR3-column mapping.
#' @field v_gene_cols Named chain-to-V-gene-column mapping.
#' @field j_gene_cols Named chain-to-J-gene-column mapping.
#' @field db_file tcrdist3 germline database.
#' @field compute_distances Logical initialization strategy for `TCRrep`.
#' @field store_all_cdr Logical CDR-matrix retention setting.
#' @field deduplicate Logical receptor deduplication setting.
#' @field cpus Number of CPUs used by tcrdist3.
#' @field radius Clustering edge radius.
#' @field python Resolved Python executable.
#' @field run_args Additional `TCRrep` constructor arguments.
Tcrdist3ToolConfig <- R6::R6Class(
  classname = "Tcrdist3ToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    organism = NULL,
    chains = NULL,
    cdr3_cols = NULL,
    v_gene_cols = NULL,
    j_gene_cols = NULL,
    db_file = NULL,
    compute_distances = NULL,
    store_all_cdr = NULL,
    deduplicate = NULL,
    cpus = NULL,
    radius = NULL,
    python = NULL,
    run_args = NULL,
    #' @description
    #' Create a validated tcrdist3 configuration. Users should call
    #' [conf_tcrdist3()] instead.
    #' @param values Named list of configuration values.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
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
      checkmate::assert_flag(values$compute_distances)
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


#' Run tcrdist3 analyses
#'
#' `Tcrdist3Tool` connects a [Tcrdist3ToolConfig] to input data and the Python
#' tcrdist3 runner. Create it with [tool_tcrdist3()], then call `$run()` for
#' each input. The tool itself contains neither the input nor results; every
#' execution returns a separate [Tcrdist3Run].
#'
#' See [tcrdist3] for the complete workflow, input formats, single- versus
#' paired-chain configuration, and clustering algorithm.
#'
#' @section Read-only fields:
#' \describe{
#'   \item{`name`}{Always `"tcrdist3"`.}
#'   \item{`config`}{The read-only [Tcrdist3ToolConfig] used by default for
#'   runs.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{`$help()`}{Print, but do not execute, the uv commands needed to
#'   create the managed Python environment and install tcrdist3.}
#'   \item{`$check()`}{Check that the configured Python exists and can import
#'   `tcrdist`. Report the Python path and version. Return `TRUE` when the check
#'   succeeds and `FALSE` otherwise.}
#'   \item{`$run(input, config = NULL, output_dir = NULL)`}{Normalize a text
#'   file or [immundata::ImmunData], execute tcrdist3, cluster the combined
#'   distance matrix, and return a [Tcrdist3Run]. Pass `config` to override the
#'   tool's default configuration for this run only. If `output_dir` is
#'   omitted, file inputs use the user home and persisted `ImmunData` inputs
#'   use the current snapshot's artifact directory.}
#'   \item{`$citation()`}{Reserved for citation metadata; currently reports
#'   that citation support is not implemented.}
#' }
#'
#' @section Runs and results:
#' `$run()` creates a new run directory and returns an immutable
#' [Tcrdist3Run] with:
#' \describe{
#'   \item{`status`}{`"success"` after the external process completes.}
#'   \item{`tool`}{The tool name, `"tcrdist3"`.}
#'   \item{`path`}{The run directory.}
#'   \item{`config`}{A plain-list snapshot of the configuration used.}
#'   \item{`input`}{A provenance-safe summary of the input.}
#'   \item{`artifacts`}{Paths to the generated files described below.}
#' }
#'
#' Call `run$results()` to obtain immutable [Tcrdist3Results]. Its
#' `$clusters`, `$mapping`, and `$edges` fields are lazy duckplyr tables, and
#' `$distances` is the NumPy matrix path. Call `results$annotate(idata, ...)`
#' to join `tcrdist_cluster_id` and `tcrdist_cluster_size` to `ImmunData`.
#' For an `ImmunData` run, result `input_id` is the original
#' `imd_receptor_id`, so one receptor-level assignment is propagated to all
#' annotation rows belonging to that receptor. Additional arguments are
#' passed to [immundata::annotate_receptors()].
#'
#' @section Artifacts:
#' `run$artifacts` contains:
#' \describe{
#'   \item{`input`}{The normalized table passed to the runner. It contains
#'   `input_id` and the tcrdist3 chain columns.}
#'   \item{`config`}{The resolved runner configuration as JSON.}
#'   \item{`clusters`}{One row per original TCR with `input_id`, normalized
#'   chain columns, `tcrdist_cluster_id`, and `tcrdist_cluster_size`.}
#'   \item{`mapping`}{Mapping from each original `input_id` to the
#'   deduplicated receptor used for distance calculation. Several original IDs
#'   can map to one deduplicated receptor.}
#'   \item{`edges`}{Pairs of TCRs whose combined distance is at most
#'   `config$radius`, expanded to original input IDs.}
#'   \item{`distances`}{Combined dense pairwise distance matrix in original
#'   input order, in NumPy `.npy` format.}
#'   \item{`stdout`, `stderr`}{Captured Python process logs.}
#' }
#'
#' @param config A [Tcrdist3ToolConfig] returned by [conf_tcrdist3()]. This
#'   becomes the tool's default configuration; `$run(config = ...)` can replace
#'   it for one execution.
#'
#' @inherit tcrdist3 examples
#'
#' @seealso [tcrdist3], [conf_tcrdist3()], [tool_tcrdist3()]
#'
Tcrdist3Tool <- R6::R6Class(
  classname = "Tcrdist3Tool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a tcrdist3 tool. Users should call [tool_tcrdist3()] instead.
    #' @param config A [Tcrdist3ToolConfig] returned by [conf_tcrdist3()].
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

      unix_uv <- "curl -LsSf https://astral.sh/uv/install.sh | sh"
      unix_env <- paste0(
        'uv venv "',
        "${IMMUNDAPTER_DATA:-$HOME/immundapter_data}",
        '/environments/tcrdist3"'
      )
      unix_tcrdist3 <- paste0(
        'uv pip install --python "',
        "${IMMUNDAPTER_DATA:-$HOME/immundapter_data}",
        '/environments/tcrdist3/bin/python" tcrdist3==0.3'
      )
      windows_uv <- paste0(
        'powershell -ExecutionPolicy ByPass -c "',
        'irm https://astral.sh/uv/install.ps1 | iex"'
      )
      windows_data <- paste0(
        "$data = if ($env:IMMUNDAPTER_DATA) ",
        "{ $env:IMMUNDAPTER_DATA } else ",
        '{ "$HOME\\immundapter_data" }'
      )
      windows_env <- paste0(
        "uv venv ",
        '"$data\\environments\\tcrdist3"'
      )
      windows_tcrdist3 <- paste0(
        "uv pip install --python ",
        '"$data\\environments\\tcrdist3\\Scripts\\python.exe" ',
        "tcrdist3==0.3"
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
      platforms <- cli::cli_ul(.close = FALSE)

      cli::cli_li("{.strong macOS and Linux}")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Install {.pkg uv}.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {unix_uv}}")
      cli::cli_end(commands)
      cli::cli_li("Create the Python environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {unix_env}}")
      cli::cli_end(commands)
      cli::cli_li("Install {.pkg tcrdist3} into the environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {unix_tcrdist3}}")
      cli::cli_end(commands)
      cli::cli_end(steps)

      cli::cli_li("{.strong Windows PowerShell}")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Install {.pkg uv}.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {windows_uv}}")
      cli::cli_end(commands)
      cli::cli_li("Set the data-directory variable for this session.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {windows_data}}")
      cli::cli_end(commands)
      cli::cli_li("Create the Python environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {windows_env}}")
      cli::cli_end(commands)
      cli::cli_li("Install {.pkg tcrdist3} into the environment.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {windows_tcrdist3}}")
      cli::cli_end(commands)
      cli::cli_end(steps)

      cli::cli_end(platforms)
      cli::cli_end(div)

      cli::cli_inform(c(
        "i" = paste0(
          "Set {.envvar IMMUNDAPTER_DATA} before running the commands ",
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
          organism = config$organism,
          chains = config$chains,
          db_file = config$db_file,
          compute_distances = config$compute_distances,
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

      list(
        input = input_path,
        config = config_path,
        clusters = clusters_path,
        mapping = mapping_path,
        edges = edges_path,
        distances = distances_path,
        stdout = stdout_path,
        stderr = stderr_path
      )
    }
  )
)
