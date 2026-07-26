#' Configure or create a RAxML-NG adapter
#'
#' `conf_raxml_ng()` validates RAxML-NG parameters. `tool_raxml_ng()` creates
#' the adapter from that configuration.
#'
#' @param sequence_col Column containing aligned nucleotide or amino-acid
#'   sequences.
#' @param family_col Column containing B-cell family identifiers.
#' @param tip_col Column containing stable identifiers for phylogenetic tree
#'   tips.
#' @param locus B-cell receptor locus to analyze.
#' @param analysis Tree-building RAxML-NG analysis mode: `"search"`, `"all"`,
#'   or `"fast"`.
#' @param model Evolutionary model passed to `--model`.
#' @param starting_tree Optional starting tree specification passed to `--tree`.
#' @param bootstrap_replicates Optional number of bootstrap trees used by
#'   `"all"` analysis.
#' @param threads Number of threads.
#' @param seed Optional random seed.
#' @param executable RAxML-NG executable name or path.
#' @param run_args Additional named RAxML-NG arguments.
#' @param config A [RaxmlNgToolConfig] created by `conf_raxml_ng()`.
#'
#' @return `conf_raxml_ng()` returns a RAxML-NG configuration.
#'   `tool_raxml_ng()` returns a RAxML-NG adapter skeleton.
#'
#' @seealso \url{https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model}
#' @export
conf_raxml_ng <- function(sequence_col = "sequence_alignment",
                          family_col = "imd_clone_id",
                          tip_col = immundata::imd_schema("receptor"),
                          locus = "IGH",
                          analysis = "search",
                          model = "GTR+G",
                          starting_tree = NULL,
                          bootstrap_replicates = NULL,
                          threads = 1L,
                          seed = NULL,
                          executable = "raxml-ng",
                          run_args = list()) {
  RaxmlNgToolConfig$new(list(
    sequence_col = sequence_col,
    family_col = family_col,
    tip_col = tip_col,
    locus = locus,
    analysis = analysis,
    model = model,
    starting_tree = starting_tree,
    bootstrap_replicates = bootstrap_replicates,
    threads = threads,
    seed = seed,
    executable = executable,
    run_args = run_args
  ))
}


#' @rdname conf_raxml_ng
#' @export
tool_raxml_ng <- function(config = conf_raxml_ng()) {
  RaxmlNgTool$new(config)
}


#' RAxML-NG tool configuration
#'
#' @keywords internal
RaxmlNgToolConfig <- R6::R6Class(
  classname = "RaxmlNgToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    sequence_col = NULL,
    family_col = NULL,
    tip_col = NULL,
    locus = NULL,
    analysis = NULL,
    model = NULL,
    starting_tree = NULL,
    bootstrap_replicates = NULL,
    threads = NULL,
    seed = NULL,
    executable = NULL,
    run_args = NULL,
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_string(values$sequence_col)
      checkmate::assert_string(values$family_col)
      checkmate::assert_string(values$tip_col)
      checkmate::assert_choice(values$locus, c("IGH", "IGK", "IGL"))
      checkmate::assert_choice(
        values$analysis,
        c("search", "all", "fast")
      )
      checkmate::assert_string(values$model)
      checkmate::assert_string(values$starting_tree, null.ok = TRUE)
      checkmate::assert_int(
        values$bootstrap_replicates,
        lower = 1L,
        null.ok = TRUE
      )
      checkmate::assert_int(values$threads, lower = 1L)
      checkmate::assert_int(values$seed, null.ok = TRUE)
      checkmate::assert_string(values$executable)
      checkmate::assert_list(values$run_args, names = "unique")

      super$initialize(values)
    }
  )
)


#' RAxML-NG adapter skeleton
#'
#' @keywords internal
RaxmlNgTool <- R6::R6Class(
  classname = "RaxmlNgTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    initialize = function(config = conf_raxml_ng()) {
      checkmate::assert_r6(config, classes = "RaxmlNgToolConfig")
      super$initialize("raxml-ng", config)
    }
  ),
  private = list(
    help_impl = function() {
      releases <- "https://codeberg.org/amkozlov/raxml-ng/releases"
      source <- "https://github.com/amkozlov/raxml-ng"
      brew <- "brew install raxml-ng"
      wsl <- "wsl --install"

      cli::cli_inform(c("i" = "Install RAxML-NG:"))
      div <- cli::cli_div(theme = list(
        ol = list("margin-left" = 2),
        ul = list("margin-left" = 2)
      ))
      platforms <- cli::cli_ul(.close = FALSE)

      cli::cli_li("{.strong macOS}")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Choose an installation method.")
      details <- cli::cli_ul(.close = FALSE)
      cli::cli_li("Homebrew: {.code {brew}}")
      cli::cli_li("Download the macOS binary: {.url {releases}}")
      cli::cli_end(details)
      cli::cli_end(steps)

      cli::cli_li("{.strong Linux}")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Download the Linux binary.")
      details <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.url {releases}}")
      cli::cli_end(details)
      cli::cli_end(steps)

      cli::cli_li("{.strong Windows}")
      steps <- cli::cli_ol(.close = FALSE)
      cli::cli_li("Install Windows Subsystem for Linux.")
      details <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {wsl}}")
      cli::cli_end(details)
      cli::cli_li("Use the Linux binary inside WSL.")
      details <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.url {releases}}")
      cli::cli_end(details)
      cli::cli_end(steps)

      cli::cli_end(platforms)
      cli::cli_end(div)

      cli::cli_inform(c(
        "i" = "Source builds and complete instructions:"
      ))
      div <- cli::cli_div(theme = list(
        ul = list("margin-left" = 2)
      ))
      links <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.url {source}}")
      cli::cli_end(links)
      cli::cli_end(div)
    },
    check_impl = function() {
      path <- Sys.which(private$.config$executable)
      installed <- nzchar(path)

      if (installed) {
        version <- system2(
          path,
          "--version",
          stdout = TRUE,
          stderr = TRUE
        )
        version <- version[nzchar(version)][[1L]]
        cli::cli_inform(c(
          "v" = "{version}",
          ">" = "Executable: {.path {path}}"
        ))
      } else {
        executable <- private$.config$executable
        cli::cli_inform(c(
          "x" = "RAxML-NG executable {.code {executable}} was not found.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      installed
    },
    citation_impl = function() {
      rlang::abort("RAxML-NG `citation()` is not implemented yet.")
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
      rlang::abort("RAxML-NG `run()` is not implemented yet.")
    },
    results_impl = function(run) {
      rlang::abort("RAxML-NG `results()` is not implemented yet.")
    },
    annotate_impl = function(idata, results, ...) {
      rlang::abort("RAxML-NG `annotate()` is not implemented yet.")
    }
  )
)
