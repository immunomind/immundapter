#' Calculate T-cell receptor distances and clusters with tcrdist3
#'
#' tcrdist3 measures similarity between T-cell receptors from their CDR3
#' amino-acid sequences and V/J genes. Similar receptors have smaller distance
#' values. In cluster mode, receptors connected within the selected distance
#' radius are placed in the same group.
#'
#' @section Basic workflow:
#' A tcrdist3 analysis has four steps:
#'
#' 1. Use [conf_tcrdist3()] to select the organism, receptor chains, input
#'    columns, and analysis mode.
#' 2. Use [tool_tcrdist3()] to create the tool.
#' 3. Use `run <- tool$run(input)` to analyze the data.
#' 4. Use `results <- run$results()` to read distances and clusters.
#'
#' In short:
#'
#' ```
#' config -> tool -> run -> results
#'                              |
#'                              +-> annotated ImmunData (cluster mode)
#' ```
#'
#' @section Choosing single- or paired-chain distance:
#' `chains` controls which receptor chains contribute to the distance:
#'
#' - **`chains = "beta"`:** Use only the beta chain. This also works with paired
#'   TRA/TRB `ImmunData`; the beta chain is selected from each receptor.
#' - **`chains = c("alpha", "beta")`:** Calculate alpha and beta distances and
#'   add them. The input must contain both chains.
#'
#' The tool checks that the selected chains are present in `ImmunData`.
#'
#' @section Input to `tool$run()`:
#' `tool$run()` accepts either an [immundata::ImmunData] dataset or one text
#' file. In both cases, each receptor must provide a CDR3 amino-acid sequence,
#' a V gene, and a J gene for every chain selected in `chains`.
#'
#' For `ImmunData`, the tool prepares one row per distinct
#' receptor. `cdr3_cols`, `v_gene_cols`, and `j_gene_cols` give the names of the
#' receptor features. For paired TRA/TRB data, both chains can use the same
#' column names because `ImmunData` already records which chain is TRA or TRB:
#'
#' ```
#' cdr3_cols  = c(alpha = "junction_aa", beta = "junction_aa")
#' v_gene_cols = c(alpha = "v_call", beta = "v_call")
#' j_gene_cols = c(alpha = "j_call", beta = "j_call")
#' ```
#'
#' A text input may be CSV, TSV, or tab-delimited TXT. Each row must represent
#' one receptor. A paired-chain file must have one row per paired receptor and
#' separate CDR3, V-gene, and J-gene columns for each chain. If the
#' file contains an `imd_receptor_id` column, that column is used as the
#' receptor ID. Otherwise, `cell_id` is used. If neither column is present,
#' row numbers are used.
#'
#' @section What happens during a run:
#' The tool performs these steps:
#'
#' 1. Rename the selected columns to the names required by tcrdist3.
#' 2. If `deduplicate = TRUE`, calculate identical receptors once.
#' 3. Calculate a distance matrix for each selected chain. For multi-chain
#'    analysis, add the matrices.
#' 4. Return distances in the same order as the original receptors.
#' 5. In cluster mode, connect receptor pairs with distance less than or equal
#'    to `radius`. Each connected group becomes one cluster.
#'
#' A group can be connected through intermediate receptors. Therefore, two
#' receptors in the same group can have a direct distance greater than
#' `radius` when other receptors form a path between them.
#'
#' @section Output from `tool$run()`:
#' `tool$run()` returns a [Tcrdist3Run]. It records where the results were
#' saved. To read them, run `results <- run$results()`.
#'
#' Both analysis modes provide:
#'
#' - **`results$mapping`:** Shows which original receptors were identical and
#'   therefore calculated once. Use `dplyr::collect()` to read the table into R.
#' - **`results$distances`:** Path to the pairwise distance matrix. The file is
#'   in Python NumPy `.npy` format; it is not an R table.
#'
#' Cluster mode also provides:
#'
#' - **`results$clusters`:** One row per original receptor, with its cluster
#'   number and cluster size.
#' - **`results$edges`:** Receptor pairs with distance less than or equal to
#'   `radius`.
#' - **`results$annotate(idata)`:** A new `ImmunData` dataset with
#'   `tcrdist_cluster_id` and `tcrdist_cluster_size` receptor columns.
#'
#' Distance mode does not create clusters and cannot annotate `ImmunData`.
#'
#' @section Managed Python environment:
#' When `python = NULL`, the tool uses its own Python environment under
#' `~/immundapter_data/environments/tcrdist3`. On Windows, it uses the
#' corresponding `Scripts/python.exe` path. Set `IMMUNDAPTER_DATA` before
#' creating the settings to use another root directory. `tool$help()`
#' prints the commands that create this environment and install the supported
#' tcrdist3 version. These commands are instructions; `$help()` does not run
#' them.
#'
#' @param mode Result type. `"cluster"` creates receptor groups and allows the
#'   original `ImmunData` to be annotated. `"distance"` creates only the
#'   distance matrix and the table linking it to the input receptors.
#' @param organism Host organism: `"human"` or `"mouse"`. It must match the V/J
#'   gene names in the input.
#' @param chains One or two unique chains selected from `"alpha"`, `"beta"`,
#'   `"gamma"`, and `"delta"`. These are passed to `TCRrep` and determine which
#'   chain distances contribute to the combined distance matrix.
#' @param cdr3_cols Named character vector that gives the CDR3 amino-acid column
#'   for every selected chain, for example `c(beta = "junction_aa")`.
#' @param v_gene_cols Named character vector that gives the V-gene column for
#'   every selected chain. Gene names and alleles must be recognized by the
#'   selected tcrdist3 germline database.
#' @param j_gene_cols Named character vector that gives the J-gene column for
#'   every selected chain.
#' @param db_file Germline database passed to
#'   `tcrdist.repertoire.TCRrep`. tcrdist3 includes the default database.
#' @param store_all_cdr Whether tcrdist3 keeps separate distance matrices for
#'   individual CDR regions in addition to its per-chain matrices.
#' @param deduplicate Whether identical receptors are calculated once. Output
#'   is still reported for every original receptor.
#' @param cpus Positive integer number of CPU cores used by tcrdist3.
#' @param radius Largest combined distance at which two receptors are directly
#'   connected during clustering. It must be a non-negative integer and is used
#'   only when `mode = "cluster"`.
#' @param python Python executable name or path. `NULL` selects the managed
#'   tcrdist3 environment described above.
#' @param run_args Additional named arguments passed to `TCRrep`. This is for
#'   advanced tcrdist3 settings. Values must be valid JSON values. Do not
#'   include arguments controlled by immundapter, such as `cell_df`, `chains`,
#'   or `index_cols`.
#'
#' @param config A [Tcrdist3ToolConfig] returned by [conf_tcrdist3()].
#'
#' @return [conf_tcrdist3()] returns a read-only
#'   [Tcrdist3ToolConfig]. [tool_tcrdist3()] returns a reusable
#'   [Tcrdist3Tool].
#'
#' @examplesIf interactive()
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
#' # The same settings can run on ImmunData.
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
#' # 5. Join receptor-level `tcrdist_cluster_id` assignments back.
#' annotated <- idata_results$annotate(idata)
#'
#' # Distance-only runs create the input links and distance matrix.
#' distance_tool <- tool_tcrdist3(conf_tcrdist3(mode = "distance"))
#' distance_run <- distance_tool$run(
#'   idata,
#'   output_dir = file.path(tempdir(), "tcrdist3-distance")
#' )
#' distance_results <- distance_run$results()
#' distance_results$distances
#' # distance_results$annotate(idata) gives an actionable error because there
#' # are no receptor-level annotations in a distance-only result.
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
#' @name tcrdist3
#' @seealso [Tcrdist3Tool], [immundata::ImmunData],
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


#' @rdname tcrdist3
#' @export
tool_tcrdist3 <- function(config = conf_tcrdist3()) {
  Tcrdist3Tool$new(config)
}


#' @rdname tcrdist3
#'
#' @description
#' `Tcrdist3ToolConfig` stores validated tcrdist3 settings. It contains no
#' receptor data and no results. Create it with [conf_tcrdist3()] and pass it
#' to [tool_tcrdist3()]. The settings cannot be changed after creation; call
#' [conf_tcrdist3()] again when different settings are needed.
#'
#' @field mode Whether the run produces clusters or distances only.
#' @field organism Organism used for V-gene and J-gene interpretation.
#' @field chains TCR chains included in the distance.
#' @field cdr3_cols CDR3 amino-acid column for each selected chain.
#' @field v_gene_cols V-gene column for each selected chain.
#' @field j_gene_cols J-gene column for each selected chain.
#' @field db_file tcrdist3 germline database.
#' @field store_all_cdr Whether to keep separate CDR distance matrices.
#' @field deduplicate Whether to calculate identical receptors once.
#' @field cpus Number of CPUs used by tcrdist3.
#' @field radius Clustering edge radius.
#' @field python Python program used to run tcrdist3.
#' @field run_args Additional advanced settings passed to `TCRrep`.
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


#' Run tcrdist3 analyses
#'
#' `Tcrdist3Tool` connects a [Tcrdist3ToolConfig] to input data and the Python
#' tcrdist3 runner. Create it with [tool_tcrdist3()], then call `$run()` for
#' each input. The tool itself contains neither the input nor results; every
#' run returns a separate [Tcrdist3Run].
#'
#' See [tcrdist3] for the complete workflow, input formats, single- versus
#' paired-chain settings, and clustering method.
#'
#' @section Read-only fields:
#'
#' - **`name`:** Always `"tcrdist3"`.
#' - **`config`:** The [Tcrdist3ToolConfig] used by default for new runs.
#'
#' @section Methods:
#'
#' - **`$help()`:** Print the commands needed to install tcrdist3. The commands
#'   are not run automatically.
#' - **`$check()`:** Check whether Python and tcrdist3 are ready. Return `TRUE`
#'   when they are ready and `FALSE` otherwise.
#' - **`$run(input, config = NULL, output_dir = NULL)`:** Analyze a text file or
#'   [immundata::ImmunData] dataset and return a [Tcrdist3Run]. Pass `config` to
#'   use different settings for this run only.
#' - **`$citation()`:** Reserved for future citation support.
#'
#' @section Runs and results:
#' `$run()` creates a new output directory and returns a [Tcrdist3Run]. Its
#' `$path` gives the output directory, `$config` gives the settings used, and
#' `$artifacts` gives the paths to files created by tcrdist3.
#'
#' Call `run$results()` to access the results. Both
#' [Tcrdist3ClusterResults] and [Tcrdist3DistanceResults] provide `$mapping`
#' and the `$distances` NumPy matrix path. Cluster results also provide
#' `$clusters` and `$edges` tables. Their
#' `$annotate(idata, ...)` joins `tcrdist_cluster_id` and
#' `tcrdist_cluster_size`.
#' For an `ImmunData` run, result `input_id` is the original
#' `imd_receptor_id`, so one cluster is added to every row for that receptor.
#' Additional arguments are passed
#' to [immundata::annotate_receptors()].
#'
#' Distance results cannot add clusters to `ImmunData`. Run the analysis with
#' `mode = "cluster"` before using `$annotate()`.
#'
#' @section Files created:
#' `run$artifacts` contains:
#'
#' - **`input`:** Table sent to tcrdist3.
#' - **`config`:** Settings sent to Python, stored as JSON.
#' - **`clusters`:** Cluster mode only. One row per original TCR, with its
#'   cluster number and cluster size.
#' - **`mapping`:** Shows which identical receptors were calculated once.
#' - **`edges`:** Cluster mode only. Receptor pairs with distance less than or
#'   equal to `config$radius`.
#' - **`distances`:** Pairwise distance matrix in NumPy `.npy` format.
#' - **`stdout`, `stderr`:** Python messages and errors.
#'
#' @param config A [Tcrdist3ToolConfig] returned by [conf_tcrdist3()]. This
#'   becomes the tool's default settings; `$run(config = ...)` can replace
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
