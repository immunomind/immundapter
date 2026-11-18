#' Infer family trees with RAxML-NG
#'
#' RAxML-NG builds a phylogenetic tree for each sequence family. This
#' tool takes aligned nucleotide sequences from [immundata::ImmunData] and
#' runs each family independently. It does not define families, align sequences,
#' or reconstruct germlines. Complete those steps before running RAxML-NG.
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
#' `additional_tips` optionally supplies backend-only tips that are appended
#' after the observed tips are extracted. It must be a data frame or duckplyr
#' table with `family_id`, `tip_id`, and `aligned_sequence` columns. Its family
#' values must already occur in the extracted input. In a partitioned analysis,
#' an optional `partition_spec` column must agree with the observed family
#' specification; when omitted, that specification is inherited. Additional
#' tips are never added to or joined into `ImmunData`.
#'
#' `outgroups` optionally identifies one or more existing tips per family. It
#' must be a data frame or duckplyr table with `family_id` and `tip_id` columns.
#' The identified tips may come from `ImmunData` or `additional_tips`. RAxML-NG
#' output is rooted by the adapter, which writes a rooted copy of each
#' successful tree while preserving the raw RAxML-NG tree.
#'
#' `locus` selects the receptor chain used for the tree. With paired
#' receptors, only the selected chain is used. The selected locus must agree
#' with the chains stored in the `ImmunData` dataset.
#'
#' If `partition_col` is set, it must contain one consistent partition
#' specification per family in the form `name=start-end`, with partitions
#' separated by semicolons. For example, `IGH=1-300;IGK=301-540` covers a
#' 540-site alignment. Partitions must be unique, non-overlapping, and cover
#' every alignment site.
#'
#' @section Analysis modes:
#' The tool supports two RAxML-NG tree-inference commands:
#'
#' - **`"search"`:** Standard maximum-likelihood tree search.
#' - **`"fast"`:** Faster tree search.
#'
#' Set `ancestral = TRUE` to run RAxML-NG ancestral-state reconstruction after
#' each successful family tree. This produces a node-labelled tree, reconstructed
#' internal-node sequences, and marginal state probabilities. A failed
#' ancestral stage does not discard the successful family tree. Other modes
#' that use existing trees or produce other types of output, including
#' `evaluate`, `bootstrap`, and `support`, are not part of this tool. The
#' default nucleotide model is `GTR+G`. Set `seed` and `threads` for
#' reproducible execution, or omit them to retain RAxML-NG's defaults.
#' Partitioned analyses default to scaled branch-length linkage.
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
#' - **`results$failed_families`:** Families that failed, with their exit status
#'   and diagnostic message.
#' - **`results$tip_map`:** The original value from `tip_col` and the shorter
#'   tip label used in the RAxML-NG files.
#' - **`results$outgroups`:** The original and internal outgroup tip names.
#' - **`results$ancestral_index`:** Paths to successful ancestral-state files.
#' - **`results$failed_ancestral`:** Families whose optional ancestral stage
#'   failed, with diagnostics.
#' - **`results$read_trees()`:** One Newick tree per family, returned as a named
#'   character vector. When outgroups were supplied, these are rooted copies.
#' - **`results$read_ancestral_trees()`:** Node-labelled ancestral trees.
#' - **`results$read_ancestral_states()`:** A duckplyr table of reconstructed
#'   internal-node sequences.
#'
#' Tree files remain in the run directory. Trees cannot be added as columns to
#' `ImmunData`. Use `$read_trees()`, `$tree_index`, or `$tip_map` for downstream
#' phylogenetic analysis.
#'
#' Each family has its own directory containing the prepared FASTA file,
#' optional partition file, RAxML-NG output, and log files. A failed family does
#' not discard successful results. The returned run has status `"partial"` when
#' one or more tree or requested ancestral stages fail.
#'
#' @param alignment_col Name of the column containing the aligned
#'   nucleotide sequence for each tree tip. This is not IgBLAST's AIRR
#'   `sequence_alignment` field, which aligns each sequence to a germline rather
#'   than to the other sequences in its family.
#' @param family_col Column containing B-cell clone, lineage, or
#'   family names or numbers. Each value is analyzed independently.
#' @param tip_col Name of the column containing the sequence name used at each
#'   tree tip. The default is the `ImmunData` receptor ID column.
#' @param locus Receptor locus to use. It must agree with the receptor chains
#'   declared by the `ImmunData` schema.
#' @param analysis Tree-inference mode. `"search"` performs the standard
#'   maximum-likelihood search; `"fast"` performs the faster search.
#' @param model Non-empty RAxML-NG nucleotide model string. The same model is
#'   used for every partition in a partitioned analysis.
#' @param seed Optional positive integer random seed. `NULL` retains the
#'   RAxML-NG default.
#' @param threads Optional positive integer thread count. `NULL` retains the
#'   RAxML-NG automatic default.
#' @param ancestral Whether to run RAxML-NG ancestral-state reconstruction on
#'   every successfully inferred family tree.
#' @param partition_col Optional column containing one partition specification
#'   per family.
#' @param brlen Optional branch-length linkage mode: `"linked"`, `"scaled"`,
#'   or `"unlinked"`. This setting requires `partition_col`; partitioned runs
#'   default to `"scaled"`.
#' @param config A [RaxmlNgToolConfig] created by [conf_raxml_ng()].
#'
#' @return `conf_raxml_ng()` returns a read-only [RaxmlNgToolConfig].
#'   `tool_raxml_ng()` returns a reusable [RaxmlNgTool].
#'
#' @examplesIf identical(Sys.getenv("IMD_RUN_EXAMPLES"), "true")
#' sequences <- c(
#'   "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
#'   "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
#'   "ACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGT",
#'   "ACGTACGTACGTACGTACGTACGTACGTACGTTCGTACGT"
#' )
#' input_file <- tempfile(fileext = ".tsv")
#' utils::write.table(data.frame(
#'   cell_id = paste0("cell-", seq_along(sequences)),
#'   locus = "IGH", sequence = sequences, clone_id = "family-a",
#'   clone_sequence_alignment = sequences, duplicate_count = 1L
#' ), input_file, sep = "\t", row.names = FALSE, quote = FALSE)
#' idata <- immundata::read_repertoires(
#'   path = input_file,
#'   schema = immundata::make_receptor_schema(
#'     features = "sequence", chains = "IGH"
#'   ),
#'   barcode_col = "cell_id", locus_col = "locus",
#'   umi_col = "duplicate_count", repertoire_schema = NULL,
#'   preprocess = NULL, postprocess = NULL, rename_columns = NULL,
#'   output_folder = tempfile("raxml-immundata-")
#' )
#' config <- conf_raxml_ng(
#'   alignment_col = "clone_sequence_alignment",
#'   family_col = "clone_id",
#'   analysis = "search"
#' )
#' tool <- tool_raxml_ng(config)
#' tool$check()
#'
#' run <- tool$run(idata, output_dir = tempfile("raxml-ng-run-"))
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
                          analysis = c("search", "fast"),
                          model = "GTR+G",
                          seed = NULL,
                          threads = NULL,
                          ancestral = FALSE,
                          partition_col = NULL,
                          brlen = NULL) {
  analysis <- match.arg(analysis)
  if (!is.null(partition_col) && is.null(brlen)) {
    brlen <- "scaled"
  }

  RaxmlNgToolConfig$new(list(
    alignment_col = alignment_col,
    family_col = family_col,
    tip_col = tip_col,
    locus = locus,
    analysis = analysis,
    model = model,
    seed = seed,
    threads = threads,
    ancestral = ancestral,
    partition_col = partition_col,
    brlen = brlen
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
#' @field locus Receptor locus analyzed.
#' @field analysis RAxML-NG inference mode.
#' @field model RAxML-NG nucleotide model.
#' @field seed Optional RAxML-NG random seed.
#' @field threads Optional RAxML-NG thread count.
#' @field ancestral Whether to reconstruct internal-node ancestral states.
#' @field partition_col Optional partition-specification column.
#' @field brlen Optional partition branch-length linkage mode.
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
    model = NULL,
    seed = NULL,
    threads = NULL,
    ancestral = NULL,
    partition_col = NULL,
    brlen = NULL,
    #' @description
    #' Create checked RAxML-NG settings. Users should call
    #' [conf_raxml_ng()] instead.
    #' @param values A named list of RAxML-NG settings.
    initialize = function(values) {
      checkmate::assert_string(values$alignment_col, min.chars = 1L)
      checkmate::assert_string(values$family_col, min.chars = 1L)
      checkmate::assert_string(values$tip_col, min.chars = 1L)
      checkmate::assert_string(
        values$partition_col,
        min.chars = 1L,
        null.ok = TRUE
      )
      if (anyDuplicated(c(
        values$alignment_col,
        values$family_col,
        values$tip_col,
        values$partition_col
      ))) {
        rlang::abort(
          paste(
            "`alignment_col`, `family_col`, `tip_col`, and `partition_col`",
            "must be different."
          )
        )
      }
      checkmate::assert_string(values$locus, min.chars = 1L)
      checkmate::assert_choice(values$analysis, c("search", "fast"))
      checkmate::assert_string(values$model, min.chars = 1L)
      if (grepl("[[:space:]]", values$model)) {
        rlang::abort("`model` cannot contain whitespace.")
      }
      checkmate::assert_count(values$seed, positive = TRUE, null.ok = TRUE)
      checkmate::assert_count(values$threads, positive = TRUE, null.ok = TRUE)
      checkmate::assert_flag(values$ancestral)
      checkmate::assert_choice(
        values$brlen,
        c("linked", "scaled", "unlinked"),
        null.ok = TRUE
      )
      if (is.null(values$partition_col) && !is.null(values$brlen)) {
        rlang::abort("`brlen` requires `partition_col`.")
      }

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
#' - **`$run(input, config = NULL, output_dir = NULL, additional_tips = NULL,
#'   outgroups = NULL)`:**
#'   Prepare one FASTA alignment per family, run RAxML-NG, and return a
#'   [RaxmlNgRun].
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
    },
    #' @description
    #' Run RAxML-NG with optional backend-only tips.
    #' @param input An [immundata::ImmunData] dataset.
    #' @param config Optional [RaxmlNgToolConfig] overriding the tool default.
    #' @param output_dir Optional output directory.
    #' @param additional_tips Optional data frame or duckplyr table with
    #'   `family_id`, `tip_id`, `aligned_sequence`, and optionally
    #'   `partition_spec` columns.
    #' @param outgroups Optional data frame or duckplyr table with `family_id`
    #'   and `tip_id` columns identifying tips on which to root each family
    #'   tree.
    run = function(input,
                   config = NULL,
                   output_dir = NULL,
                   additional_tips = NULL,
                   outgroups = NULL) {
      super$run(
        input,
        config = config,
        output_dir = output_dir,
        additional_tips = additional_tips,
        outgroups = outgroups
      )
    }
  ),
  private = list(
    help_impl = function() {
      releases <- "https://codeberg.org/amkozlov/raxml-ng/releases"
      source <- "https://github.com/amkozlov/raxml-ng"
      cli::cli_inform(c(
        "i" = "Install RAxML-NG:",
        "*" = "macOS: {.code brew install raxml-ng}, or download the binary from {.url {releases}}.",
        "*" = "Linux: download the binary from {.url {releases}}.",
        "*" = "Windows: run {.code wsl --install}, then use the Linux binary from {.url {releases}}.",
        "i" = "Source builds and complete instructions: {.url {source}}"
      ))
    },
    check_impl = function() {
      path <- raxml_ng_executable()
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
    validate_input_impl = function(input,
                                   config,
                                   additional_tips = NULL,
                                   outgroups = NULL) {
      if (!inherits(input, "ImmunData")) {
        rlang::abort(
          "RAxML-NG input must be an `immundata::ImmunData` object."
        )
      }
      input
    },
    run_impl = function(input,
                        config,
                        run_dir,
                        additional_tips = NULL,
                        outgroups = NULL) {
      executable <- raxml_ng_executable()
      if (!nzchar(executable)) {
        rlang::abort("The RAxML-NG executable `raxml-ng` was not found.")
      }

      prepared <- raxml_ng_collect_input(input, config, additional_tips)
      outgroups <- raxml_ng_prepare_outgroups(outgroups, prepared)
      input_path <- fs::path(run_dir, "raxml-ng-input.tsv")
      config_path <- fs::path(run_dir, "raxml-ng-config.json")
      families_path <- fs::path(run_dir, "raxml-ng-families.tsv")
      tree_index_path <- fs::path(run_dir, "raxml-ng-tree-index.tsv")
      failed_path <- fs::path(run_dir, "raxml-ng-failed-families.tsv")
      tip_map_path <- fs::path(run_dir, "raxml-ng-tip-map.tsv")
      outgroups_path <- fs::path(run_dir, "raxml-ng-outgroups.tsv")
      ancestral_index_path <- fs::path(
        run_dir,
        "raxml-ng-ancestral-index.tsv"
      )
      failed_ancestral_path <- fs::path(
        run_dir,
        "raxml-ng-failed-ancestral.tsv"
      )
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
      outgroup_maps <- vector("list", length(family_ids))
      normalized <- vector("list", length(family_ids))

      for (family_index in seq_along(family_ids)) {
        family_id <- family_ids[[family_index]]
        family <- prepared[prepared$family_id == family_id, , drop = FALSE]
        family$raxml_tip <- sprintf("tip_%06d", seq_len(nrow(family)))
        family_outgroups <- outgroups[
          outgroups$family_id == family_id,
          ,
          drop = FALSE
        ]
        outgroup_tips <- family$raxml_tip[
          match(family_outgroups$tip_id, family$tip_id)
        ]
        family_name <- sprintf("family-%06d", family_index)
        family_dir <- fs::path(families_dir, family_name)
        fs::dir_create(family_dir)
        alignment_path <- fs::path(family_dir, "alignment.fasta")
        model_path <- if (is.null(config$partition_col)) {
          NULL
        } else {
          fs::path(family_dir, "partitions.txt")
        }
        prefix <- fs::path(family_dir, "result")
        stdout_path <- fs::path(family_dir, "stdout.log")
        stderr_path <- fs::path(family_dir, "stderr.log")
        best_tree <- paste0(prefix, ".raxml.bestTree")
        rooted_tree <- fs::path(family_dir, "result.rootedTree")
        log_path <- paste0(prefix, ".raxml.log")
        alignment_lengths <- unique(nchar(family$aligned_sequence))

        tip_maps[[family_index]] <- data.frame(
          family_id = family_id,
          raxml_tip = family$raxml_tip,
          tip_id = family$tip_id,
          is_outgroup = family$raxml_tip %in% outgroup_tips,
          stringsAsFactors = FALSE
        )
        outgroup_maps[[family_index]] <- data.frame(
          family_id = rep(family_id, length(outgroup_tips)),
          tip_id = family_outgroups$tip_id,
          raxml_tip = outgroup_tips,
          stringsAsFactors = FALSE
        )
        normalized[[family_index]] <- family[c(
          "family_id", "tip_id", "raxml_tip", "aligned_sequence"
        )]

        fasta <- as.vector(rbind(
          paste0(">", family$raxml_tip),
          family$aligned_sequence
        ))
        writeLines(fasta, alignment_path, useBytes = TRUE)
        if (!is.null(model_path)) {
          partitions <- raxml_ng_parse_partitions(
            family$partition_spec[[1L]],
            alignment_lengths[[1L]]
          )
          writeLines(
            sprintf(
              "%s, %s = %.0f-%.0f",
              config$model,
              partitions$name,
              partitions$start,
              partitions$end
            ),
            model_path,
            useBytes = TRUE
          )
        }
        args <- raxml_ng_command_args(
          config,
          alignment_path,
          prefix,
          model_path
        )

        command <- paste(c(shQuote(executable), args), collapse = " ")
        exit_status <- tryCatch(
          raxml_ng_execute(
            executable,
            args,
            stdout_path,
            stderr_path
          ),
          error = function(condition) {
            writeLines(conditionMessage(condition), stderr_path)
            1L
          }
        )

        status <- "failed"
        error <- NA_character_
        if (isTRUE(exit_status == 0L) && fs::file_exists(best_tree)) {
          status <- "success"
        } else {
          error <- raxml_ng_execution_error(
            exit_status,
            c(stderr_path, stdout_path),
            "RAxML-NG did not produce a best-tree file."
          )
        }
        if (status == "success" && length(outgroup_tips) > 0L) {
          root_error <- tryCatch(
            {
              raxml_ng_root_tree(best_tree, rooted_tree, outgroup_tips)
              NULL
            },
            error = function(condition) conditionMessage(condition)
          )
          if (!is.null(root_error)) {
            status <- "failed"
            error <- gsub(
              "[\t\r\n]+",
              " ",
              paste("Tree rooting failed:", root_error)
            )
          }
        }

        tree_path <- if (status == "success") {
          if (length(outgroup_tips) > 0L) rooted_tree else best_tree
        } else {
          NA_character_
        }

        ancestral_status <- if (config$ancestral) {
          "skipped"
        } else {
          "not_requested"
        }
        ancestral_prefix <- fs::path(family_dir, "ancestral")
        ancestral_stdout <- fs::path(family_dir, "ancestral-stdout.log")
        ancestral_stderr <- fs::path(family_dir, "ancestral-stderr.log")
        raw_ancestral_tree <- paste0(
          ancestral_prefix,
          ".raxml.ancestralTree"
        )
        rooted_ancestral_tree <- fs::path(
          family_dir,
          "ancestral.rootedTree"
        )
        ancestral_states <- paste0(
          ancestral_prefix,
          ".raxml.ancestralStates"
        )
        ancestral_probabilities <- paste0(
          ancestral_prefix,
          ".raxml.ancestralProbs"
        )
        ancestral_log <- paste0(ancestral_prefix, ".raxml.log")
        ancestral_command <- NA_character_
        ancestral_exit_status <- NA_integer_
        ancestral_error <- NA_character_
        ancestral_tree <- NA_character_

        if (config$ancestral && status == "success") {
          ancestral_args <- raxml_ng_command_args(
            config,
            alignment_path,
            ancestral_prefix,
            model_path,
            best_tree
          )
          ancestral_command <- paste(
            c(shQuote(executable), ancestral_args),
            collapse = " "
          )
          ancestral_exit_status <- tryCatch(
            raxml_ng_execute(
              executable,
              ancestral_args,
              ancestral_stdout,
              ancestral_stderr
            ),
            error = function(condition) {
              writeLines(conditionMessage(condition), ancestral_stderr)
              1L
            }
          )
          expected <- c(
            raw_ancestral_tree,
            ancestral_states,
            ancestral_probabilities
          )
          if (isTRUE(ancestral_exit_status == 0L) &&
              all(fs::file_exists(expected))) {
            ancestral_status <- "success"
            ancestral_tree <- raw_ancestral_tree
          } else {
            ancestral_status <- "failed"
            ancestral_error <- raxml_ng_execution_error(
              ancestral_exit_status,
              c(ancestral_stderr, ancestral_stdout),
              "RAxML-NG did not produce all ancestral-state files."
            )
          }
          if (ancestral_status == "success" && length(outgroup_tips) > 0L) {
            root_error <- tryCatch(
              {
                raxml_ng_root_tree(
                  raw_ancestral_tree,
                  rooted_ancestral_tree,
                  outgroup_tips
                )
                NULL
              },
              error = function(condition) conditionMessage(condition)
            )
            if (is.null(root_error)) {
              ancestral_tree <- rooted_ancestral_tree
            } else {
              ancestral_status <- "failed"
              ancestral_tree <- NA_character_
              ancestral_error <- gsub(
                "[\t\r\n]+",
                " ",
                paste("Ancestral-tree rooting failed:", root_error)
              )
            }
          }
        }

        family_records[[family_index]] <- data.frame(
          family_id = family_id,
          family_index = family_index,
          status = status,
          n_tips = nrow(family),
          n_outgroups = length(outgroup_tips),
          alignment_length = alignment_lengths,
          alignment = as.character(alignment_path),
          model_file = if (is.null(model_path)) {
            NA_character_
          } else {
            as.character(model_path)
          },
          best_tree = raxml_ng_existing_path(best_tree),
          rooted_tree = raxml_ng_existing_path(rooted_tree),
          tree = tree_path,
          raxml_log = raxml_ng_existing_path(log_path),
          stdout = as.character(stdout_path),
          stderr = as.character(stderr_path),
          command = command,
          exit_status = as.integer(exit_status),
          error = error,
          ancestral_status = ancestral_status,
          ancestral_tree = ancestral_tree,
          raw_ancestral_tree = raxml_ng_existing_path(raw_ancestral_tree),
          ancestral_states = raxml_ng_existing_path(ancestral_states),
          ancestral_probabilities = raxml_ng_existing_path(
            ancestral_probabilities
          ),
          ancestral_log = raxml_ng_existing_path(ancestral_log),
          ancestral_stdout = if (config$ancestral && status == "success") {
            as.character(ancestral_stdout)
          } else {
            NA_character_
          },
          ancestral_stderr = if (config$ancestral && status == "success") {
            as.character(ancestral_stderr)
          } else {
            NA_character_
          },
          ancestral_command = ancestral_command,
          ancestral_exit_status = as.integer(ancestral_exit_status),
          ancestral_error = ancestral_error,
          stringsAsFactors = FALSE
        )
      }

      normalized <- do.call(rbind, normalized)
      family_records <- do.call(rbind, family_records)
      tip_maps <- do.call(rbind, tip_maps)
      outgroup_maps <- do.call(rbind, outgroup_maps)
      tree_index <- family_records[
        family_records$status == "success",
        c(
          "family_id", "family_index", "n_tips", "n_outgroups",
          "alignment_length", "tree", "best_tree", "rooted_tree",
          "raxml_log"
        ),
        drop = FALSE
      ]
      failed <- family_records[
        family_records$status == "failed",
        c(
          "family_id", "family_index", "n_tips", "alignment_length",
          "alignment", "model_file", "best_tree", "rooted_tree", "stdout",
          "stderr", "command", "exit_status", "error"
        ),
        drop = FALSE
      ]
      ancestral_index <- family_records[
        family_records$ancestral_status == "success",
        c(
          "family_id", "family_index", "ancestral_tree",
          "raw_ancestral_tree", "ancestral_states",
          "ancestral_probabilities", "ancestral_log"
        ),
        drop = FALSE
      ]
      failed_ancestral <- family_records[
        family_records$ancestral_status == "failed",
        c(
          "family_id", "family_index", "alignment", "best_tree",
          "ancestral_stdout", "ancestral_stderr", "ancestral_command",
          "ancestral_exit_status", "ancestral_error"
        ),
        drop = FALSE
      ]

      utils::write.table(normalized, file = input_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(family_records, file = families_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(tree_index, file = tree_index_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(failed, file = failed_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(tip_maps, file = tip_map_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(outgroup_maps, file = outgroups_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(ancestral_index, file = ancestral_index_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(failed_ancestral, file = failed_ancestral_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

      artifacts <- list(
        input = input_path,
        config = config_path,
        families = families_path,
        tree_index = tree_index_path,
        failed_families = failed_path,
        tip_map = tip_map_path,
        outgroups = outgroups_path,
        ancestral_index = ancestral_index_path,
        failed_ancestral = failed_ancestral_path
      )
      attr(artifacts, "run_status") <- if (
        nrow(failed) > 0L || nrow(failed_ancestral) > 0L
      ) {
        "partial"
      } else {
        "success"
      }
      artifacts
    }
  )
)
