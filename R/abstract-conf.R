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
    #' @description
    #' Create a configuration from validated values. Concrete configuration
    #' classes validate their values before they call this method.
    #' @param values A named list of configuration values. Every name must
    #'   match a public field of the concrete configuration class.
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

    #' @description
    #' Return a plain-list copy of the configuration.
    #' @return A named list containing the configuration fields and their
    #'   current values.
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
