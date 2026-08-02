#' Infer B-cell family trees with RAxML-NG
#'
#' The RAxML-NG adapter reads family-level aligned sequences from
#' [immundata::ImmunData] annotations, writes one FASTA alignment per family,
#' and independently infers one maximum-likelihood tree per family.
#' It does not assign clones, reconstruct germlines, or align sequences.
#'
#' @section Input contract:
#' `tool$run()` accepts an `ImmunData` object. `alignment_col`, `family_col`,
#' and `tip_col` name columns in `input$annotations`. They are annotations, not
#' receptor-schema features: a family alignment is derived from a group of
#' receptors and must not define receptor identity.
#'
#' `family_col` groups receptors into independent alignments. `tip_col`
#' identifies the leaves within each family, and `alignment_col` contains one
#' family-aligned nucleotide string for each leaf. Repeated annotation rows
#' with the same family, tip, and alignment are collapsed. Conflicting rows for
#' one family and tip are rejected. All entries in a family must have the same
#' aligned length.
#'
#' For a paired receptor schema, `locus` filters the long annotation table to
#' the selected chain. For a single-chain schema, the configured locus must
#' match the schema chain when the schema declares one.
#'
#' @section Analysis modes:
#' The initial adapter supports two RAxML-NG tree-inference commands:
#' \describe{
#'   \item{`"search"`}{A standard maximum-likelihood tree search.}
#'   \item{`"fast"`}{The RAxML-NG fast-search shortcut.}
#' }
#'
#' Modes that consume existing trees or produce different scientific objects,
#' including `evaluate`, `bootstrap`, `support`, and `ancestral`, are not part
#' of this interface. Both supported modes use the nucleotide model `GTR+G`.
#' RAxML-NG controls its own threads, starting trees, random seed, and other
#' execution settings.
#'
#' @section Runs and results:
#' Every family is isolated in its own directory with its normalized FASTA,
#' RAxML-NG files, stdout, and stderr. The run succeeds only when every family
#' produces its expected best tree. If a family fails, all completed family
#' outputs and logs remain in the run directory and `run()` reports an error.
#'
#' `run$results()` returns [RaxmlNgResults]. Its `families`, `tree_index`, and
#' `tip_map` fields are lazy duckplyr tables. Tree files remain on disk; use
#' `results$read_trees()` to read selected Newick strings.
#'
#' @param alignment_col Annotation column containing the nucleotide sequence
#'   aligned within its family. This is not IgBLAST's AIRR
#'   `sequence_alignment` field.
#' @param family_col Annotation column containing B-cell clone, lineage, or
#'   family identifiers. Each value is analyzed independently.
#' @param tip_col Stable identifier for tree tips. The default is the canonical
#'   `ImmunData` receptor identifier.
#' @param locus B-cell receptor locus to analyze.
#' @param analysis Tree-inference mode: `"search"` or `"fast"`.
#' @param config A [RaxmlNgToolConfig] created by [conf_raxml_ng()].
#'
#' @return `conf_raxml_ng()` returns a read-only [RaxmlNgToolConfig].
#'   `tool_raxml_ng()` returns a reusable [RaxmlNgTool].
#'
#' @examples
#' \dontrun{
#' config <- conf_raxml_ng(
#'   alignment_col = "clone_sequence_alignment",
#'   family_col = "clone_id",
#'   analysis = "search"
#' )
#' tool <- tool_raxml_ng(config)
#' tool$check()
#'
#' run <- tool$run(idata, output_dir = "analysis/raxml-ng")
#' results <- run$results()
#' dplyr::collect(results$tree_index)
#' trees <- results$read_trees()
#' }
#'
#' @seealso
#'   \url{https://github.com/amkozlov/raxml-ng/wiki/Input-data},
#'   \url{https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings}
#' @export
conf_raxml_ng <- function(alignment_col = "clone_sequence_alignment",
                          family_col = "clone_id",
                          tip_col = immundata::imd_schema("receptor"),
                          locus = "IGH",
                          analysis = c("search", "fast")) {
  analysis <- match.arg(analysis)

  RaxmlNgToolConfig$new(list(
    alignment_col = alignment_col,
    family_col = family_col,
    tip_col = tip_col,
    locus = locus,
    analysis = analysis
  ))
}


#' @rdname conf_raxml_ng
#' @export
tool_raxml_ng <- function(config = conf_raxml_ng()) {
  RaxmlNgTool$new(config)
}


#' RAxML-NG tool configuration
#'
#' `RaxmlNgToolConfig` is an immutable recipe for one family-level RAxML-NG
#' analysis. Create it with [conf_raxml_ng()].
#'
#' @field alignment_col Family-aligned nucleotide annotation column.
#' @field family_col Family identifier annotation column.
#' @field tip_col Stable tree-tip identifier.
#' @field locus B-cell receptor locus analyzed.
#' @field analysis RAxML-NG inference mode.
#'
#' @keywords internal
RaxmlNgToolConfig <- R6::R6Class(
  classname = "RaxmlNgToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    alignment_col = NULL,
    family_col = NULL,
    tip_col = NULL,
    locus = NULL,
    analysis = NULL,
    #' @description
    #' Create a validated RAxML-NG configuration. Users should call
    #' [conf_raxml_ng()] instead.
    #' @param values A named list of RAxML-NG configuration values.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_string(values$alignment_col, min.chars = 1L)
      checkmate::assert_string(values$family_col, min.chars = 1L)
      checkmate::assert_string(values$tip_col, min.chars = 1L)
      if (anyDuplicated(c(
        values$alignment_col,
        values$family_col,
        values$tip_col
      ))) {
        rlang::abort(
          "`alignment_col`, `family_col`, and `tip_col` must be different."
        )
      }
      checkmate::assert_choice(values$locus, c("IGH", "IGK", "IGL"))
      checkmate::assert_choice(values$analysis, c("search", "fast"))

      super$initialize(values)
    }
  )
)


#' Run RAxML-NG family tree inference
#'
#' `RaxmlNgTool` connects a [RaxmlNgToolConfig] to an `ImmunData` input and
#' returns a separate [RaxmlNgRun] for every execution.
#'
#' @section Methods:
#' \describe{
#'   \item{`$help()`}{Print platform-specific installation guidance.}
#'   \item{`$check()`}{Check and report the installed executable and version.}
#'   \item{`$run(input, config = NULL, output_dir = NULL)`}{Prepare aligned
#'   family FASTAs, run each family, and return a [RaxmlNgRun].}
#'   \item{`$citation()`}{Reserved for citation metadata.}
#' }
#'
#' @keywords internal
RaxmlNgTool <- R6::R6Class(
  classname = "RaxmlNgTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a RAxML-NG adapter. Users should call [tool_raxml_ng()] instead.
    #' @param config A [RaxmlNgToolConfig] returned by [conf_raxml_ng()].
    initialize = function(config = conf_raxml_ng()) {
      checkmate::assert_r6(config, classes = "RaxmlNgToolConfig")
      super$initialize("raxml-ng", config, RaxmlNgRun)
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
      path <- Sys.which("raxml-ng")
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
        cli::cli_inform(c(
          "x" = "RAxML-NG executable {.code raxml-ng} was not found.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      installed
    },
    citation_impl = function() {
      rlang::abort("RAxML-NG `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      if (!inherits(input, "ImmunData")) {
        rlang::abort(
          "RAxML-NG input must be an `immundata::ImmunData` object."
        )
      }
      checkmate::assert_r6(input, classes = "ImmunData")
      input
    },
    run_impl = function(input, config, run_dir) {
      executable <- Sys.which("raxml-ng")
      if (!nzchar(executable)) {
        rlang::abort("The RAxML-NG executable `raxml-ng` was not found.")
      }

      prepared <- raxml_ng_collect_input(input, config)
      input_path <- fs::path(run_dir, "raxml-ng-input.tsv")
      config_path <- fs::path(run_dir, "raxml-ng-config.json")
      families_path <- fs::path(run_dir, "raxml-ng-families.tsv")
      tree_index_path <- fs::path(run_dir, "raxml-ng-tree-index.tsv")
      failed_path <- fs::path(run_dir, "raxml-ng-failed-families.tsv")
      tip_map_path <- fs::path(run_dir, "raxml-ng-tip-map.tsv")
      families_dir <- fs::path(run_dir, "families")
      fs::dir_create(families_dir)

      jsonlite::write_json(
        config$as_list(),
        path = config_path,
        auto_unbox = TRUE,
        null = "null",
        pretty = TRUE
      )

      family_ids <- sort(unique(prepared$family_id))
      family_records <- vector("list", length(family_ids))
      tip_maps <- vector("list", length(family_ids))
      normalized <- vector("list", length(family_ids))

      for (family_index in seq_along(family_ids)) {
        family_id <- family_ids[[family_index]]
        family <- prepared[prepared$family_id == family_id, , drop = FALSE]
        family$raxml_tip <- sprintf("tip_%06d", seq_len(nrow(family)))
        family_name <- sprintf("family-%06d", family_index)
        family_dir <- fs::path(families_dir, family_name)
        fs::dir_create(family_dir)
        alignment_path <- fs::path(family_dir, "alignment.fasta")
        prefix <- fs::path(family_dir, "result")
        stdout_path <- fs::path(family_dir, "stdout.log")
        stderr_path <- fs::path(family_dir, "stderr.log")
        best_tree <- paste0(prefix, ".raxml.bestTree")
        log_path <- paste0(prefix, ".raxml.log")
        alignment_lengths <- unique(nchar(family$aligned_sequence))

        tip_maps[[family_index]] <- data.frame(
          family_id = family_id,
          raxml_tip = family$raxml_tip,
          tip_id = family$tip_id,
          stringsAsFactors = FALSE
        )
        normalized[[family_index]] <- family[c(
          "family_id", "tip_id", "raxml_tip", "aligned_sequence"
        )]

        fasta <- as.vector(rbind(
          paste0(">", family$raxml_tip),
          family$aligned_sequence
        ))
        writeLines(fasta, alignment_path)
        args <- c(
          paste0("--", config$analysis),
          "--msa", shQuote(alignment_path),
          "--model", "GTR+G",
          "--prefix", shQuote(prefix)
        )

        command <- paste(c(executable, args), collapse = " ")
        exit_status <- tryCatch(
          suppressWarnings(system2(
            executable,
            args,
            stdout = stdout_path,
            stderr = stderr_path
          )),
          error = function(condition) {
            writeLines(conditionMessage(condition), stderr_path)
            1L
          }
        )

        status <- "failed"
        error <- NA_character_
        if (identical(exit_status, 0L) && fs::file_exists(best_tree)) {
          status <- "success"
        } else {
          log_lines <- character()
          for (path in c(stderr_path, stdout_path)) {
            if (fs::file_exists(path)) {
              log_lines <- c(log_lines, readLines(path, warn = FALSE))
            }
          }
          log_lines <- log_lines[nzchar(log_lines)]
          error <- if (length(log_lines) == 0L) {
            "RAxML-NG exited without an error message."
          } else {
            paste(tail(log_lines, 20L), collapse = " | ")
          }
        }

        family_records[[family_index]] <- data.frame(
          family_id = family_id,
          family_index = family_index,
          status = status,
          n_tips = nrow(family),
          alignment_length = alignment_lengths,
          alignment = as.character(alignment_path),
          best_tree = if (fs::file_exists(best_tree)) {
            as.character(best_tree)
          } else {
            NA_character_
          },
          raxml_log = if (fs::file_exists(log_path)) {
            as.character(log_path)
          } else {
            NA_character_
          },
          stdout = as.character(stdout_path),
          stderr = as.character(stderr_path),
          command = command,
          error = error,
          stringsAsFactors = FALSE
        )
      }

      normalized <- do.call(rbind, normalized)
      family_records <- do.call(rbind, family_records)
      tip_maps <- do.call(rbind, tip_maps)
      tree_index <- family_records[
        family_records$status == "success",
        c(
          "family_id", "family_index", "n_tips", "alignment_length",
          "best_tree", "raxml_log"
        ),
        drop = FALSE
      ]
      failed <- family_records[
        family_records$status == "failed",
        c(
          "family_id", "family_index", "n_tips", "alignment_length",
          "alignment", "stdout", "stderr", "error"
        ),
        drop = FALSE
      ]

      utils::write.table(normalized, file = input_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(family_records, file = families_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(tree_index, file = tree_index_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(failed, file = failed_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(tip_maps, file = tip_map_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

      if (nrow(failed) > 0L) {
        details <- paste(
          sprintf("%s: %s", failed$family_id, failed$error),
          collapse = "\n"
        )
        rlang::abort(paste(
          "RAxML-NG failed for one or more families.",
          sprintf("Run files were retained at: %s", run_dir),
          details,
          sep = "\n"
        ))
      }

      artifacts <- list(
        input = input_path,
        config = config_path,
        families = families_path,
        tree_index = tree_index_path,
        tip_map = tip_map_path
      )
      artifacts
    }
  )
)
