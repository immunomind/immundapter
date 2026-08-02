#' Parsed RAxML-NG results
#'
#' @field families Lazy duckplyr table with one row per requested family.
#' @field tree_index Lazy duckplyr table indexing successful best-tree files.
#' @field tip_map Lazy duckplyr mapping from RAxML-safe tip labels to original
#'   tip identifiers.
#'
#' @keywords internal
RaxmlNgResults <- R6::R6Class(
  classname = "RaxmlNgResults",
  inherit = ToolResults,
  cloneable = FALSE,
  public = list(
    #' @description Read successful trees as Newick strings.
    #' @param families Optional character vector of family identifiers. `NULL`
    #'   reads every successful tree.
    #' @return A named character vector of Newick strings.
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
          "RAxML-NG trees cannot be stored as tabular `ImmunData`",
          "annotations. Use `$tree_index`, `$tip_map`, or `$read_trees()`",
          "to work with the file-backed phylogenetic results."
        ),
        class = "immundapter_unannotatable_results"
      )
    }
  )
)
