#' Base class for tool configuration
#'
#' Tool configurations are small, validated collections of named values.
#' Concrete classes declare their fields and validate them before calling this
#' constructor.
#'
#' @keywords internal
AbstractConfig <- R6::R6Class(
  classname = "AbstractConfig",
  cloneable = FALSE,
  public = list(
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")

      unknown <- setdiff(names(values), names(self))
      if (length(unknown) > 0L) {
        rlang::abort(sprintf(
          "Unknown configuration field(s): %s.",
          paste(unknown, collapse = ", ")
        ))
      }

      for (field in names(values)) {
        self[[field]] <- values[[field]]
      }
      private$.fields <- names(values)

      for (field in private$.fields) {
        lockBinding(field, self)
      }
    },

    as_list = function() {
      stats::setNames(
        lapply(private$.fields, function(field) self[[field]]),
        private$.fields
      )
    }
  ),
  private = list(
    .fields = character()
  )
)
