#' Infer B-cell family trees with RAxML-NG
#'
#' RAxML-NG builds a phylogenetic tree for each B-cell clone or family. This
#' tool takes aligned nucleotide sequences from [immundata::ImmunData] and
#' runs each family independently. It does not identify clones, reconstruct
#' germline sequences, or align sequences. Complete these steps before running
#' RAxML-NG.
#'
#' @section Input to `tool$run()`:
#' `tool$run()` accepts one [immundata::ImmunData] dataset. It does not accept a
#' FASTA file directly. The dataset must contain three columns:
#'
#' - **`family_col`:** Clone, lineage, or family name. A separate tree is built
#'   for each value in this column.
#' - **`tip_col`:** Sequence name shown at the tips, or leaves, of the tree.
#'   Each sequence within a family must have a different name.
#' - **`alignment_col`:** Aligned nucleotide sequence. All sequences in one
#'   family must have the same length.
#'
#' All three values must be present and non-empty. Repeated rows with the same
#' family, tip, and aligned sequence are used once. A tip cannot have two
#' different sequences in the same family.
#'
#' `locus` selects the B-cell receptor chain used for the tree. With paired
#' receptors, only the selected chain is used. The selected locus must agree
#' with the chains stored in the `ImmunData` dataset.
#'
#' @section Analysis modes:
#' The tool supports two RAxML-NG tree-inference commands:
#'
#' - **`"search"`:** Standard maximum-likelihood tree search.
#' - **`"fast"`:** Faster tree search.
#'
#' Modes that use existing trees or produce other types of output,
#' including `evaluate`, `bootstrap`, `support`, and `ancestral`, are not part
#' of this tool. Both supported modes use the nucleotide model `GTR+G`.
#' RAxML-NG chooses its own number of computer threads, starting trees, random
#' seed, and other run settings.
#'
#' @section Output from `tool$run()`:
#' `tool$run()` returns a [RaxmlNgRun]. It records where the results were
#' saved. To read the trees, first run `results <- run$results()`.
#'
#' [RaxmlNgResults] provides:
#'
#' - **`results$families`:** One row per family, showing whether the analysis
#'   succeeded and where its files were saved. Use `dplyr::collect()` to read
#'   the table into R.
#' - **`results$tree_index`:** The families that produced a tree and the path
#'   to each tree file.
#' - **`results$tip_map`:** The original value from `tip_col` and the shorter
#'   tip label used in the RAxML-NG files.
#' - **`results$read_trees()`:** One Newick tree per family, returned as a named
#'   character vector. The vector names are the family values.
#'
#' Tree files remain in the run directory. Trees cannot be added as columns to
#' `ImmunData`. Use `$read_trees()`, `$tree_index`, or `$tip_map` for downstream
#' phylogenetic analysis.
#'
#' Each family has its own directory containing the prepared FASTA file,
#' RAxML-NG output, and log files. A run is successful only when every family
#' produces a best tree. If one family fails, completed files are kept and
#' `tool$run()` reports an error.
#'
#' @param alignment_col Name of the column containing the aligned
#'   nucleotide sequence for each tree tip. This is not IgBLAST's AIRR
#'   `sequence_alignment` field, which aligns each sequence to a germline rather
#'   than to the other sequences in its family.
#' @param family_col Column containing B-cell clone, lineage, or
#'   family names or numbers. Each value is analyzed independently.
#' @param tip_col Name of the column containing the sequence name used at each
#'   tree tip. The default is the `ImmunData` receptor ID column.
#' @param locus B-cell receptor locus to use: `"IGH"`, `"IGK"`, or `"IGL"`.
#' @param analysis Tree-inference mode. `"search"` performs the standard
#'   maximum-likelihood search; `"fast"` performs the faster search.
#' @param config A [RaxmlNgToolConfig] created by [conf_raxml_ng()].
#'
#' @return `conf_raxml_ng()` returns a read-only [RaxmlNgToolConfig].
#'   `tool_raxml_ng()` returns a reusable [RaxmlNgTool].
#'
#' @examplesIf interactive()
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
#' @name raxml-ng
#' @seealso
#'   [RAxML-NG input](https://github.com/amkozlov/raxml-ng/wiki/Input-data),
#'   [RAxML-NG output](https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings)
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


#' @rdname raxml-ng
#' @export
tool_raxml_ng <- function(config = conf_raxml_ng()) {
  RaxmlNgTool$new(config)
}


#' @rdname raxml-ng
#'
#' @description
#' `RaxmlNgToolConfig` stores validated settings for one type of family-level
#' RAxML-NG analysis. It contains no sequences and no trees. Create it with
#' [conf_raxml_ng()] and pass it to [tool_raxml_ng()]. The settings cannot be
#' changed after creation; call [conf_raxml_ng()] again to use different
#' settings.
#'
#' @field alignment_col Column containing family-aligned nucleotides.
#' @field family_col Column containing the family name or number.
#' @field tip_col Column containing the sequence name used at each tree tip.
#' @field locus B-cell receptor locus analyzed.
#' @field analysis RAxML-NG inference mode.
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
    #' Create checked RAxML-NG settings. Users should call
    #' [conf_raxml_ng()] instead.
    #' @param values A named list of RAxML-NG settings.
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
#' `RaxmlNgTool` uses one [RaxmlNgToolConfig] as its default settings. Create it
#' with [tool_raxml_ng()]. See [raxml-ng] for the required `ImmunData` input and
#' available tree output.
#'
#' @section Methods:
#'
#' - **`$help()`:** Print installation instructions.
#' - **`$check()`:** Check whether RAxML-NG is installed and report its version.
#' - **`$run(input, config = NULL, output_dir = NULL)`:** Prepare one FASTA
#'   alignment per family, run RAxML-NG, and return a [RaxmlNgRun].
#' - **`$citation()`:** Reserved for future citation support.
#'
#' @keywords internal
RaxmlNgTool <- R6::R6Class(
  classname = "RaxmlNgTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a RAxML-NG tool. Users should call [tool_raxml_ng()] instead.
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
