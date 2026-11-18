#' Parsed RAxML-NG results
#'
#' `RaxmlNgResults` provides tables that describe the family runs, trees, and
#' ancestral-state files. Use `$read_trees()` and
#' `$read_ancestral_trees()` to read trees as Newick strings. The tables stay on
#' disk until you use `dplyr::collect()`.
#'
#' RAxML-NG trees cannot be added as columns to `ImmunData`.
#'
#' @field families Read-only table with one row per family, including whether
#'   the analysis succeeded and where its files were saved.
#' @field tree_index Read-only table with one row per successful tree.
#' @field failed_families Read-only table with one row per failed family.
#' @field tip_map Read-only table with the original tip names and the shorter
#'   names used in RAxML-NG files.
#' @field outgroups Read-only table mapping requested outgroups to RAxML-NG tip
#'   names.
#' @field ancestral_index Read-only table with paths to successful ancestral
#'   outputs.
#' @field failed_ancestral Read-only table with one row per failed ancestral
#'   stage.
#'
#' @keywords internal
RaxmlNgResults <- R6::R6Class(
  classname = "RaxmlNgResults",
  inherit = ToolResults,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Read successful trees as Newick strings. A Newick string is a compact
    #' text representation of a phylogenetic tree.
    #' @param families Optional character vector of family names or numbers. `NULL`
    #'   reads every successful tree.
    #' @return A named character vector of Newick strings. Vector names are the
    #'   family names or numbers.
    read_trees = function(families = NULL) {
      index <- private$family_index(
        self$tree_index,
        families,
        "No successful RAxML-NG tree is available for"
      )

      trees <- vapply(index$tree, function(path) {
        checkmate::assert_file_exists(path)
        paste(readLines(path, warn = FALSE), collapse = "")
      }, character(1L))
      names(trees) <- index$family_id
      trees
    },
    #' @description
    #' Read successful node-labelled ancestral trees as Newick strings.
    #' @param families Optional character vector of family names or numbers.
    #'   `NULL` reads every successful ancestral tree.
    #' @return A named character vector of Newick strings.
    read_ancestral_trees = function(families = NULL) {
      index <- private$family_index(
        self$ancestral_index,
        families,
        "No successful RAxML-NG ancestral tree is available for"
      )

      trees <- vapply(index$ancestral_tree, function(path) {
        checkmate::assert_file_exists(path)
        paste(readLines(path, warn = FALSE), collapse = "")
      }, character(1L))
      names(trees) <- index$family_id
      trees
    },
    #' @description
    #' Read reconstructed internal-node sequences.
    #' @param families Optional character vector of family names or numbers.
    #'   `NULL` reads every successful ancestral-state file.
    #' @return A duckplyr table with `family_id`, `node_id`, and `sequence`.
    read_ancestral_states = function(families = NULL) {
      index <- private$family_index(
        self$ancestral_index,
        families,
        "No successful RAxML-NG ancestral states are available for"
      )

      states <- lapply(seq_len(nrow(index)), function(row) {
        path <- index$ancestral_states[[row]]
        checkmate::assert_file_exists(path)
        table <- utils::read.delim(
          path,
          header = FALSE,
          col.names = c("node_id", "sequence"),
          colClasses = "character",
          quote = "",
          comment.char = "",
          stringsAsFactors = FALSE
        )
        data.frame(
          family_id = index$family_id[[row]],
          node_id = table$node_id,
          sequence = table$sequence,
          stringsAsFactors = FALSE
        )
      })
      if (length(states) == 0L) {
        states <- list(data.frame(
          family_id = character(),
          node_id = character(),
          sequence = character(),
          stringsAsFactors = FALSE
        ))
      }
      do.call(rbind, states) |>
        duckplyr::as_duckdb_tibble()
    }
  ),
  active = list(
    families = function(value) {
      private$artifact("families", value)
    },
    tree_index = function(value) {
      private$artifact("tree_index", value)
    },
    failed_families = function(value) {
      private$artifact("failed_families", value)
    },
    tip_map = function(value) {
      private$artifact("tip_map", value)
    },
    outgroups = function(value) {
      private$artifact("outgroups", value)
    },
    ancestral_index = function(value) {
      private$artifact("ancestral_index", value)
    },
    failed_ancestral = function(value) {
      private$artifact("failed_ancestral", value)
    }
  ),
  private = list(
    family_index = function(index, families, unavailable) {
      checkmate::assert_character(families, unique = TRUE, null.ok = TRUE)
      index <- index |>
        dplyr::collect() |>
        as.data.frame()
      index$family_id <- as.character(index$family_id)
      if (is.null(families)) {
        return(index)
      }

      missing <- setdiff(families, index$family_id)
      if (length(missing) > 0L) {
        rlang::abort(sprintf(
          "%s family/families: %s.",
          unavailable,
          paste(missing, collapse = ", ")
        ))
      }
      index[match(families, index$family_id), , drop = FALSE]
    },
    annotate_impl = function(idata, ...) {
      rlang::abort(
        paste(
          "RAxML-NG trees cannot be added as columns to `ImmunData`.",
          "Use `$tree_index`, `$tip_map`, or `$read_trees()`",
          "to work with the phylogenetic trees and their files."
        ),
        class = "immundapter_unannotatable_results"
      )
    }
  )
)
