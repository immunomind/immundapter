#' Parsed RAxML-NG results
#'
#' `RaxmlNgResults` provides tables that describe the family runs and tree
#' files. Use `$read_trees()` to read the trees as Newick strings. The tables
#' stay on disk until you use `dplyr::collect()`.
#'
#' RAxML-NG trees cannot be added as columns to `ImmunData`.
#'
#' @field families Read-only table with one row per family, including whether
#'   the analysis succeeded and where its files were saved.
#' @field tree_index Read-only table with one row per successful tree.
#' @field tip_map Read-only table with the original tip names and the shorter
#'   names used in RAxML-NG files.
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
      checkmate::assert_character(families, unique = TRUE, null.ok = TRUE)
      index <- self$tree_index |>
        dplyr::collect() |>
        as.data.frame()
      index$family_id <- as.character(index$family_id)

      if (!is.null(families)) {
        missing <- setdiff(families, index$family_id)
        if (length(missing) > 0L) {
          rlang::abort(sprintf(
            "No successful RAxML-NG tree is available for family/families: %s.",
            paste(missing, collapse = ", ")
          ))
        }
        index <- index[match(families, index$family_id), , drop = FALSE]
      }

      trees <- vapply(index$best_tree, function(path) {
        checkmate::assert_file_exists(path)
        paste(readLines(path, warn = FALSE), collapse = "")
      }, character(1L))
      names(trees) <- index$family_id
      trees
    }
  ),
  active = list(
    families = function(value) {
      if (!missing(value)) {
        rlang::abort("`families` is read-only.")
      }
      path <- private$.run$artifacts$families
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    },
    tree_index = function(value) {
      if (!missing(value)) {
        rlang::abort("`tree_index` is read-only.")
      }
      path <- private$.run$artifacts$tree_index
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    },
    tip_map = function(value) {
      if (!missing(value)) {
        rlang::abort("`tip_map` is read-only.")
      }
      path <- private$.run$artifacts$tip_map
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    }
  ),
  private = list(
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
